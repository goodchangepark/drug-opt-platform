"""Sanitized request trace context for attributable scientific mutations."""

from __future__ import annotations

import re
import uuid
from contextvars import ContextVar, Token
from dataclasses import dataclass

from fastapi import Request

_SAFE = re.compile(r"[^A-Za-z0-9._:@/-]+")


def _clean(value: str | None, *, default: str = "", limit: int = 120) -> str:
    text = _SAFE.sub("_", str(value or "").strip())
    return text[:limit] or default


@dataclass(frozen=True)
class RequestTrace:
    caller_type: str
    caller_name: str
    request_id: str
    workflow_id: str
    execution_id: str
    conversation_id: str
    route_action: str
    transaction_id: str


_TRACE: ContextVar[RequestTrace | None] = ContextVar("drugopt_request_trace", default=None)


def trace_from_request(request: Request) -> RequestTrace:
    caller = _clean(request.headers.get("X-DrugOPT-Caller"), default="UI_OR_API")
    caller_type = "INTEGRATION" if caller.upper() not in {"UI", "UI_OR_API", "API"} else "APPLICATION"
    request_id = _clean(request.headers.get("X-Request-ID"), default=str(uuid.uuid4()))
    workflow_id = _clean(request.headers.get("X-Workflow-ID"))
    execution_id = _clean(request.headers.get("X-Execution-ID"))
    conversation_id = _clean(request.headers.get("X-Conversation-ID"))
    route_action = _clean(f"{request.method}:{request.url.path}", limit=200)
    return RequestTrace(
        caller_type=caller_type,
        caller_name=caller,
        request_id=request_id,
        workflow_id=workflow_id,
        execution_id=execution_id,
        conversation_id=conversation_id,
        route_action=route_action,
        transaction_id=str(uuid.uuid4()),
    )


def set_request_trace(trace: RequestTrace) -> Token:
    return _TRACE.set(trace)


def reset_request_trace(token: Token) -> None:
    _TRACE.reset(token)


def current_request_trace() -> RequestTrace:
    trace = _TRACE.get()
    if trace is not None:
        return trace
    return RequestTrace(
        caller_type="BACKGROUND_OR_DIRECT",
        caller_name="UNATTRIBUTED_NON_REQUEST",
        request_id="",
        workflow_id="",
        execution_id="",
        conversation_id="",
        route_action="DIRECT_CALL",
        transaction_id=str(uuid.uuid4()),
    )
