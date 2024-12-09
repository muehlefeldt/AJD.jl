using PyCall

py"""
def py_add(x, y):
    return x + y
"""

py"py_add"(1, 5)