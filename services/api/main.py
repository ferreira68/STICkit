from fastapi import FastAPI
from pydantic import BaseModel

app = FastAPI(title="STICkit API")


class EchoReq(BaseModel):
    text: str


class EchoResp(BaseModel):
    text: str


@app.post("/echo", response_model=EchoResp)
def echo(req: EchoReq) -> EchoResp:
    return EchoResp(text=req.text)
