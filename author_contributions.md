From a file called author_contributions.tsv, with columns initials and category, extract the initials of the authors per category. Categories can be comma separated, so should be split accordingly. 

Contribution categories are encoded by a single letter as follows:

- C	biospecimen collection, pathological analysis and clinical data collection
- G	genetic analysis and molecular biology studies
- M	drafted the manuscript
- S	statistics
- O	oversaw, coordinated

```python
import pandas as pd
# Load the author contributions data
df = pd.read_csv('author_contributions.tsv', sep='\t')
# Create multiple lines for authors with multiple categories
df = df.assign(category=df['category'].str.split(',')).explode('category')
# Group by category, aggregate initials into lists, and add dots between letters of initials
# e.g., "WDC" becomes "W.D.C."
contributions = df.groupby('category')['initials'].apply(lambda x: ['.'.join(i) + '.' for i in x]).to_dict()
# create a mapping of category codes to full descriptions
category_map = {
    'C': 'Biospecimen collection, pathological analysis and clinical data collection',
    'G': 'Genetic analysis and molecular biology studies',
    'M': 'Drafted the manuscript',
    'S': 'Statistics',
    'O': 'Oversaw and coordinated'
}
# Replace category codes with full descriptions
for code, initials in contributions.items():
    full_desc = category_map.get(code, code)
    print(f"{full_desc}: {', '.join(initials)},")
```