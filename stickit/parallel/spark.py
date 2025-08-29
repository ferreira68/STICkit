from __future__ import annotations
from typing import Optional

try:
    from pyspark.sql import SparkSession, DataFrame  # type: ignore
except Exception:  # pragma: no cover
    SparkSession = None  # type: ignore
    DataFrame = None  # type: ignore


def get_spark(app_name: str = "stickit") -> Optional["SparkSession"]:
    if SparkSession is None:
        return None
    return SparkSession.builder.appName(app_name).getOrCreate()


def etl_example(spark: "SparkSession", input_path: str, limit: int = 0) -> "DataFrame":
    df = spark.read.parquet(input_path)
    if limit:
        df = df.limit(limit)
    return df
