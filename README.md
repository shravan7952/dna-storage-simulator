# DNA Storage Simulator

This project is developed by **Codecell.ai** to demonstrate the principles of DNA-based data storage.  
It provides a complete simulation pipeline for encoding digital files into DNA sequences and decoding them back, including error simulation and recovery.

## Features
- Convert any file into:
  - Binary (bits)
  - DNA sequence (ATCG encoding)
- Analyze DNA sequences:
  - Base composition
  - GC content
  - Homopolymer detection
- Oligo generation with:
  - Primers
  - Reed–Solomon error correction
  - Replication
- Error simulation:
  - Substitution mutations
  - Insertions and deletions (indels)
  - Oligo dropout
- DNA Fountain (outer code) implementation
- Recovery with success metrics:
  - RS decoding
  - Fountain peeling decoder
  - Success/failure statistics and plots
- FASTA export of oligos

## Installation

It is recommended to use **Anaconda** to manage dependencies.

### Step 1: Clone the repository
```bash
git clone https://github.com/shravan7952/dna-storage-simulator.git
cd dna-storage-simulator
````

### Step 2: Create a new conda environment

```bash
conda create -n biocomp python=3.10 -y
conda activate biocomp
```

### Step 3: Install dependencies

```bash
pip install -r requirements.txt
```

Dependencies include:

* streamlit
* pandas
* reedsolo

## Usage

Run the Streamlit application:

```bash
streamlit run app.py
```

This will open the simulator in your browser.

## Deployment

This project can also be deployed on **Streamlit Cloud** for public access.
Simply connect this repository to Streamlit Cloud and deploy.

## License

This project is open-source under the MIT License.

## Credits

Developed by **Codecell.ai**
