# PCR Probabilities

This plays with conditional probabilities around PCR tests.
We are going to calculate the chance of select a random person who is
infected of COVID.
After that we need the false positive and negative rates of PCR tests.
All data is based on Madrid.

## Data source:
- https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data#field-description
- https://www.thelancet.com/journals/lanres/article/PIIS2213-2600(20)30453-7/fulltext

```python
import pandas as pd
import re
```

```python
pd.options.display.max_colwidth = 120
PTRN_EVENT = re.compile("P?\(?\s*([\w']+)(\s?\|\s?([\w']+))?\s*\)?")
DECIMALS = 4
```

```python
class P:
    """Probability"""

    def __init__(self, name, value, decimals=DECIMALS):
        self.decimals = decimals
        if decimals is not None:
            value = round(value, decimals)
        self.value = value
        main, _, base = re.search(PTRN_EVENT, name).groups()
        self.name = self.get_name(main, base)
        self.main = main
        self.base = base

    @classmethod
    def get_name(cls, main, base=None):
        if base is None:
            s = f"{main}"
        else:
            s = f"{main} | {base}"

        name = f"P({s})"
        return name

    def __repr__(self):
        return f"{self.name} = {self.value}"

    def __str__(self):
        return f"{self.name} = {self.value}"

    def c(self, name=None):
        "Calculates complemenatary"
        value = 1 - self.value

        if name is not None:
            return P(name, value, decimals=self.decimals)

        not_str = "'"
        if not_str in self.main:
            not_event = self.main.replace(not_str, "")
        else:
            not_event = f"{self.main}{not_str}"

        name = self.get_name(not_event, self.base)

        return P(name, value, decimals=self.decimals)

    def bayes(self, Pevent, Pbase):
        """Applies bayes rule.

        Examples
        --------
        >>> p_ab = P("a | b", value=.5)
        >>> p_ba = p_ab.bayes(p_b, p_a)
        """
        assert Pevent.main == self.main, "Pevent unexpected"
        assert Pbase.main == self.base, "Pbase unexpected"

        value = self.value * Pbase.value / Pevent.value
        name = self.get_name(self.base, self.main)
        return P(name, value, decimals=self.decimals)

    @staticmethod
    def Pmain(Pmain_base_list, Pbase_list):
        """Probability of main"""
        value = 0
        for Pmbi, Pbi in zip(Pmain_base_list, Pbase_list):
            assert Pmbi.base == Pbi.main, "Wrong combination"
            value += Pmbi.value * Pbi.value

        return P(Pmbi.main, value)

    def round(self, decimals):
        self.value = round(self.value, decimals)
        return self
```

# Load data from github
https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data#field-description

```python
def load_data():
    prefix = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/'
    path = 'csse_covid_19_data/csse_covid_19_daily_reports'
    date = '10-31-2020'
    url_fmt = '{prefix}/{path}/{date}.csv'
    url = url_fmt.format(prefix=prefix, path=path, date=date)
    desease_info = pd.read_csv(url)
    return desease_info
```

```python
madrid = load_data().query("Province_State == 'Madrid'")
```

```python
# Define variables
incidence_rate = madrid.Incidence_Rate.item()
```

```python
prob_desease = P("D", value=incidence_rate / 100_000)
prob_no_desease = prob_desease.c()
```

## PCR INFO
> The current rate of operational false-positive swab tests in the UK is
> unknown; preliminary estimates show it could be somewhere between 0·8% and
> 4·0%

Source : https://www.thelancet.com/journals/lanres/article/PIIS2213-2600(20)30453-7/fulltext

```python
# P(pos | D') =
prob_positive_no_desease = P("pos | D'", value=(0.8/100 + 4/100) / 2)
prob_negative_no_desease = prob_positive_no_desease.c("neg | D'")
```

> Globally, most effort so far has been invested in turnaround times and
> low test sensitivity (ie, false negatives); one systematic review reported
> false-negative rates of between 2% and 33% in repeat sample testing.

Source : https://www.thelancet.com/journals/lanres/article/PIIS2213-2600(20)30453-7/fulltext

```python
# P(neg | D) =
prob_negative_desease = P("neg | D", value=(2/100 + 33/100) / 2)
prob_negative_desease
```

```python
# P(pos | D) =
prob_positive_desease = prob_negative_desease.c("pos | D")
prob_positive_desease
```

```python
# P(neg) = P(neg | D) * P(D) + P(neg | D') * P(D') =
prob_neg = P.Pmain(
    [prob_negative_desease, prob_negative_no_desease],
    [prob_desease, prob_no_desease]
)
prob_neg
```

## Conditional Probabilities


We are going to calculate the probability of be infected when the
result of pcr is negative.

```python
# P(D | neg) = P(neg | D) * P(D) / P(neg) =
prob_negative_desease.bayes(prob_negative, prob_desease)
prob_negative_desease
```

```python
# P(pos) = P(main | base) * P(D) + P(pos | D') * P(D') =
prob_positive = P.Pmain(
    [prob_positive_desease, prob_positive_no_desease],
    [prob_desease, prob_no_desease],
)
prob_positive
```

We are going to calculate the probability of not be infected when the
result of pcr is positive.

```python
# P(D' | pos) = P(pos | D') * P(D') / P(pos) =
prob_no_desease_positive = prob_positive_no_desease.bayes(prob_positive, prob_no_desease)
prob_no_desease_positive
```
