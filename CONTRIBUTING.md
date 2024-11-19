## Adding new ranking methods
To add a new ranking method, follow these steps.
First, write a function that takes a `dual_dataset` and returns a pandas dataframe.
This is the function that computes the metric.

If you need extra arguments, write a subparser, like so:
```python
def test_method(dual_dataset, a, b):
    print(a + b)
    return dual_dataset.case

custom_parser = ArgumentParser()
custom_parser.add_argument("a")
custom_parser.add_argument("b")
```

You can configure `parser` as much as you want.

When you are done, add to `RANKING_METHODS` your method, wrapped in a `RankingMethod`:
```python
RANKING_METHODS.update(
    {
        "id_of_method": RankingMethod(
            name = "Human name of method",
            exec = test_method,
            parser = custom_parser # If None, a dummy parser is added for you, with no extra args
            desc = "A long description of the method, shown in --list-methods"
        )
    }
)

```

And you're done! The extra method will show up in the list of methods.
