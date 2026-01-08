# Structure of populations.xlsx

Based on the code analysis, here's the structure needed for each sheet:

## Sheet: img1

**Purpose:** Compares haplotype carrier frequencies between aFTLD-U patients and Controls across different cohorts

**Required columns:**
- `phenotype` - Text: "aFTLD-U" or "Controls"
- `cohort` - Text: Cohort name (e.g., "Cohort 1", "Cohort 2", "Cohort 1+2")
- `hapA` - Integer: Number of hapA carriers
- `hapB` - Integer: Number of hapB carriers
- `none` - Integer: Number of non-carriers
- `sum` - Integer: Total sample size (hapA + hapB + none)

**Example structure:**
| phenotype | cohort | hapA | hapB | none | sum |
|-----------|--------|------|------|------|-----|
| aFTLD-U | Cohort 1 | | | | |
| aFTLD-U | Cohort 1+2 | | | | |
| Controls | Cohort 1 | | | | |
| Controls | Cohort 2 | | | | |

---

## Sheet: img2

**Purpose:** Shows haplotype distribution in non-aFTLD-U populations

**Required columns:**
- `cohort` - Text: Cohort/population name
- `hapA` - Integer: Number of hapA carriers
- `hapB` - Integer: Number of hapB carriers  
- `none` - Integer: Number with neither hapA nor hapB
- `none2` - Integer: Additional category (displayed as "none+hapB")
- `sum` - Integer: Total sample size

**Example structure:**
| cohort | hapA | hapB | none | none2 | sum |
|--------|------|------|------|-------|-----|
| Population 1 | | | | | |
| Population 2 | | | | | |

---

## Sheet: img3

**Purpose:** Shows repeat classification distribution across haplotypes and phenotypes

**Required columns:**
- `haplotype` - Text: "hapA", "hapB", or "hapA+B"
- `phenotype` - Text: "aFTLD-U" or "Controls"
- `classification` - Text: Repeat classification category (e.g., ">450bp/>80%CT", etc.)
- `CT repeat` - Integer: Count of CT repeat samples
- `Short CT repeat` - Integer: Count of Short CT repeat samples
- `CnT repeat` - Integer: Count of CnT repeat samples
- `12-mer repeat` - Integer: Count of 12-mer repeat samples
- `CCCTCT repeat` - Integer: Count of CCCTCT repeat samples
- `no repeat` - Integer: Count of samples with no repeat
- `sum` - Integer: Total sample size for this row

**Example structure:**
| haplotype | phenotype | classification | CT repeat | Short CT repeat | CnT repeat | 12-mer repeat | CCCTCT repeat | no repeat | sum |
|-----------|-----------|----------------|-----------|-----------------|------------|---------------|---------------|-----------|-----|
| hapA | aFTLD-U | >450bp/>80%CT | | | | | | | |
| hapA | aFTLD-U | Other category | | | | | | | |
| hapA | Controls | >450bp/>80%CT | | | | | | | |
| hapB | aFTLD-U | >450bp/>80%CT | | | | | | | |
| hapB | Controls | >450bp/>80%CT | | | | | | | |
| hapA+B | aFTLD-U | >450bp/>80%CT | | | | | | | |
| hapA+B | Controls | >450bp/>80%CT | | | | | | | |

**Note:** The classification ">450bp/>80%CT" appears in the full version but is filtered out in the simplified version of image 3.

---

## Notes:
- All count columns should contain integers (0 or positive numbers)
- The `sum` column should equal the sum of all count columns in that row
- Make sure to save the file as `populations.xlsx` in `~/local/` directory
