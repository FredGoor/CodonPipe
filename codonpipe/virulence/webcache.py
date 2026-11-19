#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Small, dependency-light HTTP layer with an on-disk cache.

Every remote resource used by the virulence-discovery pipeline goes through
:class:`WebCache`.  This buys three properties that matter for a published
pipeline:

* **Reproducibility.**  Each response is written to disk next to a sidecar
  ``.meta.json`` recording the exact URL, query parameters, HTTP status and
  UTC access time.  A run can therefore be repeated byte-for-byte later, and
  the provenance report can state when every database was queried.
* **Offline operation.**  With ``offline=True`` nothing touches the network;
  cached payloads are replayed and a missing entry raises
  :class:`OfflineCacheMiss`.  Reviewers can re-run an analysis from the cache
  folder alone.
* **Politeness.**  Requests to a given host are spaced by ``min_interval``
  seconds and retried with exponential back-off on transient failures, so the
  public APIs are not hammered.

The module deliberately uses only ``requests`` (already a CodonPipe
dependency) and the standard library.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from typing import Any, Callable, Dict, Optional, Tuple
from urllib.parse import urlparse

try:  # requests is a CodonPipe dependency, but fail softly and informatively
    import requests
except Exception as _exc:  # pragma: no cover
    requests = None
    _REQUESTS_IMPORT_ERROR = _exc
else:
    _REQUESTS_IMPORT_ERROR = None


__all__ = [
    "WebCache",
    "OfflineCacheMiss",
    "NetworkUnavailable",
    "CACHE_FORMAT_VERSION",
]

CACHE_FORMAT_VERSION = 1

#: Transient HTTP statuses that justify a retry rather than an error.
_RETRY_STATUS = (429, 500, 502, 503, 504)


class OfflineCacheMiss(RuntimeError):
    """Raised when ``offline=True`` and the requested resource is not cached."""


class NetworkUnavailable(RuntimeError):
    """Raised when a resource could not be retrieved after all retries."""


def _safe_component(text: str) -> str:
    keep = "-_."
    out = "".join(c if (c.isalnum() or c in keep) else "_" for c in str(text))
    return out[:80] or "x"


class WebCache:
    """Cached, rate-limited, retrying HTTP GET helper.

    Parameters
    ----------
    cache_dir:
        Directory holding cached payloads.  Created on demand.
    offline:
        When true, never open a socket; replay cache or raise
        :class:`OfflineCacheMiss`.
    refresh:
        When true, ignore existing cache entries and re-fetch (the new
        response replaces the cached one).
    timeout:
        Per-request timeout in seconds.
    retries:
        Number of additional attempts after the first failure.
    min_interval:
        Minimum seconds between two requests to the same host.
    user_agent:
        Value sent in the ``User-Agent`` header.  Keep it generic and
        descriptive of the tool.
    log:
        Optional ``callable(str)`` used for progress messages.
    """

    def __init__(
        self,
        cache_dir: str,
        offline: bool = False,
        refresh: bool = False,
        timeout: float = 60.0,
        retries: int = 3,
        min_interval: float = 0.34,
        user_agent: str = "CodonPipe-virulence-discovery/1.2",
        log: Optional[Callable[[str], None]] = None,
    ) -> None:
        self.cache_dir = os.path.abspath(os.path.expanduser(str(cache_dir)))
        self.offline = bool(offline)
        self.refresh = bool(refresh)
        self.timeout = float(timeout)
        self.retries = int(retries)
        self.min_interval = float(min_interval)
        self.user_agent = str(user_agent)
        self._log = log
        self._last_call: Dict[str, float] = {}
        self._session = None
        #: Access log, consumed by the provenance report.
        self.accesses: list = []
        os.makedirs(self.cache_dir, exist_ok=True)

    # ------------------------------------------------------------------ utils
    def log(self, message: str) -> None:
        if self._log is not None:
            try:
                self._log(str(message))
            except Exception:
                pass

    @property
    def session(self):
        if self._session is None:
            if requests is None:  # pragma: no cover
                raise NetworkUnavailable(
                    "The 'requests' package is required for virulence discovery "
                    f"but could not be imported: {_REQUESTS_IMPORT_ERROR}"
                )
            self._session = requests.Session()
            self._session.headers.update({"User-Agent": self.user_agent})
        return self._session

    def _throttle(self, host: str) -> None:
        last = self._last_call.get(host)
        if last is not None:
            wait = self.min_interval - (time.time() - last)
            if wait > 0:
                time.sleep(wait)
        self._last_call[host] = time.time()

    def _paths(self, url: str, params: Optional[dict], ext: str) -> Tuple[str, str]:
        canonical = url
        if params:
            items = sorted((str(k), str(v)) for k, v in params.items())
            canonical += "?" + "&".join(f"{k}={v}" for k, v in items)
        digest = hashlib.sha1(canonical.encode("utf-8")).hexdigest()[:20]
        host = _safe_component(urlparse(url).netloc or "local")
        folder = os.path.join(self.cache_dir, host)
        os.makedirs(folder, exist_ok=True)
        stem = os.path.join(folder, f"{digest}")
        return stem + ext, stem + ".meta.json"

    def _record(self, url: str, params: Optional[dict], status: Any, source: str) -> None:
        self.accesses.append(
            {
                "url": url,
                "params": dict(params or {}),
                "status": status,
                "source": source,
                "accessed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
        )

    # ------------------------------------------------------------------- core
    def get_bytes(
        self,
        url: str,
        params: Optional[dict] = None,
        ext: str = ".bin",
        accept: Optional[str] = None,
    ) -> bytes:
        """Return the response body for ``url``, from cache when possible."""
        payload_path, meta_path = self._paths(url, params, ext)

        if os.path.isfile(payload_path) and not self.refresh:
            with open(payload_path, "rb") as fh:
                data = fh.read()
            self._record(url, params, "cache", "cache")
            return data

        if self.offline:
            raise OfflineCacheMiss(
                "Offline mode is enabled and this resource is not in the cache:\n"
                f"  {url}\n  params={params}\n"
                f"  expected cache file: {payload_path}"
            )

        host = urlparse(url).netloc
        headers = {"Accept": accept} if accept else None
        last_error: Optional[BaseException] = None
        for attempt in range(self.retries + 1):
            self._throttle(host)
            try:
                resp = self.session.get(
                    url, params=params, timeout=self.timeout, headers=headers
                )
            except Exception as exc:  # network-level failure
                last_error = exc
                if attempt < self.retries:
                    time.sleep(min(2.0 ** attempt, 8.0))
                    continue
                break
            if resp.status_code in _RETRY_STATUS and attempt < self.retries:
                time.sleep(min(2.0 ** attempt, 8.0))
                continue
            if resp.status_code >= 400:
                last_error = NetworkUnavailable(
                    f"HTTP {resp.status_code} for {resp.url}"
                )
                break
            data = resp.content
            tmp = payload_path + ".part"
            with open(tmp, "wb") as fh:
                fh.write(data)
            os.replace(tmp, payload_path)
            with open(meta_path, "w", encoding="utf-8") as fh:
                json.dump(
                    {
                        "cache_format_version": CACHE_FORMAT_VERSION,
                        "url": url,
                        "resolved_url": resp.url,
                        "params": dict(params or {}),
                        "status_code": resp.status_code,
                        "content_type": resp.headers.get("Content-Type"),
                        "bytes": len(data),
                        "fetched_utc": time.strftime(
                            "%Y-%m-%dT%H:%M:%SZ", time.gmtime()
                        ),
                    },
                    fh,
                    indent=1,
                )
            self._record(url, params, resp.status_code, "network")
            return data

        self._record(url, params, "error", "network")
        raise NetworkUnavailable(
            f"Could not retrieve {url} (params={params}) after "
            f"{self.retries + 1} attempt(s): {last_error}"
        )

    def get_text(
        self,
        url: str,
        params: Optional[dict] = None,
        ext: str = ".txt",
        accept: Optional[str] = None,
    ) -> str:
        return self.get_bytes(url, params, ext=ext, accept=accept).decode(
            "utf-8", errors="replace"
        )

    def get_json(
        self, url: str, params: Optional[dict] = None, ext: str = ".json"
    ) -> Any:
        text = self.get_text(url, params, ext=ext, accept="application/json")
        return json.loads(text)

    # -------------------------------------------------------------- reporting
    def summary(self) -> Dict[str, int]:
        """Counts of cache hits versus network fetches, for the report."""
        hits = sum(1 for a in self.accesses if a["source"] == "cache")
        net = sum(1 for a in self.accesses if a["source"] == "network")
        return {"requests": len(self.accesses), "from_cache": hits, "from_network": net}
