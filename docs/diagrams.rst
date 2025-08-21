Mermaid demo
============

.. mermaid::

   flowchart TD
     A[Start] --> B{Valid?}
     B -- yes --> C[Do thing]
     B -- no  --> D[Handle error]
     C --> E[End]
