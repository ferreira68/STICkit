from __future__ import annotations
import seaborn as sns
def apply_house_style(context: str = "talk", style: str = "whitegrid") -> None:
    sns.set_theme(context=context, style=style)
