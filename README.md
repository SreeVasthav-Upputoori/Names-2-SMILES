
# 🧪 Chemical Name to SMILES Converter

A powerful, multithreaded Python tool to **convert chemical names to SMILES (Simplified Molecular Input Line Entry System)** using public chemical databases. Supports both **command-line execution** and **interactive mode**, with support for **batch processing**, **caching**, and **real-time validation**.

---

## 🚀 Features

- 🔍 Fetches SMILES from:
  - [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
  - [OPSIN](https://opsin.ch.cam.ac.uk/)
  - [NCI CACTUS](https://cactus.nci.nih.gov/)
- 📂 Input: `.csv`, `.xlsx`, or `.xls` files with chemical names
- 📄 Output: `.xlsx` file with SMILES and validity check
- ⚡ Multithreaded API requests (customizable workers and batch size)
- 💾 Caching mechanism to avoid repeated lookups
- 📈 Summary statistics and logging
- ✅ SMILES structure validation
- 🧠 Intelligent cleaning of chemical names

---

## 📦 Installation

1. Clone the repository:

```bash
git clone https://github.com/yourusername/name2smiles-converter.git
cd name2smiles-converter
```

2. Install dependencies:

```bash
pip install -r requirements.txt
```

_**Note:** Requires Python 3.7+._

---

## 🛠 Usage

### 🔹 Option 1: Interactive Mode

```bash
python main.py
```

Follow the prompts to select input file, output name, and column containing chemical names.

---

### 🔹 Option 2: Command-Line Mode

```bash
python main.py input_file.xlsx -o output_file.xlsx -c "Name"
```

| Argument         | Description                                |
|------------------|--------------------------------------------|
| `input_file`     | Input CSV or Excel file                    |
| `-o` / `--output`| Output Excel file (default: auto-named)    |
| `-c` / `--column`| Column name with chemical names (default: `Name`) |
| `--max-workers`  | Number of threads (default: 10)            |
| `--batch-size`   | Batch size per API request (default: 5)    |
| `--timeout`      | Request timeout in seconds (default: 15)   |
| `--cache-file`   | Cache file name (default: `smiles_cache.json`) |
| `--clear-cache`  | Clears saved SMILES cache                  |
| `--interactive` or `-i` | Forces interactive mode             |

---

## 🧪 Example

Input (`compounds.xlsx`):

| Name              |
|-------------------|
| Aspirin           |
| Caffeine          |
| Paracetamol       |

Output (`results.xlsx`):

| Name       | SMILES                     | SMILES_Valid |
|------------|-----------------------------|--------------|
| Aspirin    | CC(=O)OC1=CC=CC=C1C(=O)O     | True         |
| Caffeine   | CN1C=NC2=C1C(=O)N(C(=O)N2C)C | True         |
| Paracetamol| CC(=O)NC1=CC=C(C=C1)O        | True         |

---

## 📊 Output Summary

After execution, the tool prints:
- Number of compounds processed
- Number of SMILES found
- SMILES validation stats
- Cache usage
- API success/failure breakdown

---

## 🧰 Dependencies

- `pandas`
- `requests`
- `openpyxl` (for Excel I/O)
- `argparse`
- `concurrent.futures`
- `logging`

Install them with:

```bash
pip install pandas requests openpyxl
```

---

## 📁 File Structure

```
├── main.py                # Main script
├── smiles_cache.json      # Local cache (auto-created)
├── name2smiles.log        # Logs
└── requirements.txt       # Optional dependency list
```

---

## 🧠 Notes

- The tool uses multiple APIs to improve reliability.
- Caching is automatic and helps reduce redundant API calls.
- Supports thousands of molecules with efficient threading.

---

## 🤝 License

MIT License. Feel free to use, modify, and distribute.

---

## 🙋‍♀️ Author

Created by **Sree Vasthav Upputoori**  
📧 sreevasthav.upputoori@gmail.com

---

## 🌟 Acknowledgements

Thanks to the developers of:
- PubChem PUG REST API
- OPSIN (Cambridge)
- NCI CACTUS chemical resolver
