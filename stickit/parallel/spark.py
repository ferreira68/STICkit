from __future__ import annotations
from typing import Optional, TYPE_CHECKING

try:
    from pyspark.sql import SparkSession, DataFrame
except Exception as e:  # pragma: no cover
    print(f"Unable to import pyspark: {e}")
    if not TYPE_CHECKING:
        SparkSession = None
        DataFrame = None


def get_spark(app_name: str = "stickit") -> Optional["SparkSession"]:
    if SparkSession is None:
        return None
    return SparkSession.builder.appName(app_name).getOrCreate()


def etl_example(spark: "SparkSession", input_path: str, limit: int = 0) -> "DataFrame":
    df = spark.read.parquet(input_path)
    if limit:
        df = df.limit(limit)
    return df
