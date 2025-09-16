# app.py
import streamlit as st
import pandas as pd
import random
from collections import Counter, defaultdict, deque
import reedsolo
import math
import hashlib
import json
import itertools
import os

# -------------------------
# Branding (CodeCell.ai)
# -------------------------
COMPANY_NAME = "codecell.ai"
LOGO_PATH = "codecell_logo.png"  # put your logo file next to this app.py (PNG recommended)
AUTHOR = "CodeCell.ai — DNA Storage Lab"

# page config with simple icon
st.set_page_config(page_title=f"{COMPANY_NAME} — DNA Storage Simulator", layout="wide", page_icon="🧬")

# show top header with logo (if available)
col1, col2 = st.columns([1, 8])
with col1:
    if os.path.exists(LOGO_PATH):
        st.image(LOGO_PATH, width=120)
    else:
        # fallback small placeholder
        st.markdown("###")
with col2:
    st.markdown(f"## **{COMPANY_NAME}** — DNA Storage Simulator")
    st.markdown(f"**{AUTHOR}**  \nFull pipeline demo: Bits → DNA → Oligos → ECC → Fountain → Errors → Decode")

# small separator
st.write("---")

# -------------------------
# Helpers: Bits <-> Bytes
# -------------------------
def bytes_to_bits(b: bytes) -> str:
    return ''.join(f"{byte:08b}" for byte in b)

def bits_to_bytes_exact(bits: str) -> bytes:
    assert len(bits) % 8 == 0, "bit string must be multiple of 8"
    return bytes(int(bits[i:i+8], 2) for i in range(0, len(bits), 8))

# -------------------------
# Lookup mapping: 2 bits -> 2 bases (safe)
# -------------------------
LOOKUP_MAP = {'00': 'AT', '01': 'CG', '10': 'TA', '11': 'GC'}
INV_LOOKUP = {v: k for k, v in LOOKUP_MAP.items()}

def bits_to_dna_lookup(bits: str) -> str:
    if len(bits) % 2 != 0:
        bits += '0'
    return ''.join(LOOKUP_MAP[bits[i:i+2]] for i in range(0, len(bits), 2))

def dna_to_bits_lookup(dna: str) -> str:
    trim = len(dna) - (len(dna) % 2)
    dna = dna[:trim]
    return ''.join(INV_LOOKUP.get(dna[i:i+2], '00') for i in range(0, len(dna), 2))

# -------------------------
# DNA statistics helpers
# -------------------------
def base_counts(dna: str) -> dict:
    c = Counter(dna)
    return {b: c.get(b, 0) for b in ['A', 'T', 'C', 'G']}

def gc_content(dna: str) -> float:
    return 100.0 * (dna.count('G') + dna.count('C')) / len(dna) if dna else 0.0

def find_homopolymers(dna: str, min_run: int = 6) -> list:
    runs = []
    if not dna:
        return runs
    cur_base, cur_start, cur_len = dna[0], 0, 1
    for i, ch in enumerate(dna[1:], start=1):
        if ch == cur_base:
            cur_len += 1
        else:
            if cur_len >= min_run:
                runs.append((cur_base, cur_len, cur_start, i-1))
            cur_base, cur_len, cur_start = ch, 1, i
    if cur_len >= min_run:
        runs.append((cur_base, cur_len, cur_start, len(dna)-1))
    return runs

# -------------------------
# Primers & RS helper
# -------------------------
PRIMER_FWD = "ACGTACGTACGT"
PRIMER_REV = "TGCATGCATGCA"

def make_rs(nsym: int):
    return reedsolo.RSCodec(nsym)

# -------------------------
# Simple DNA Fountain (LT-code style)
# -------------------------
def split_into_blocks(payload_bytes: bytes, block_size: int):
    blocks = [payload_bytes[i:i+block_size] for i in range(0, len(payload_bytes), block_size)]
    if len(blocks) > 0 and len(blocks[-1]) < block_size:
        blocks[-1] = blocks[-1].ljust(block_size, b'\x00')
    return blocks

def fountain_generate_droplet(blocks, seed):
    rng = random.Random(seed)
    K = len(blocks)
    r = rng.random()
    if r < 0.2:
        d = 1
    elif r < 0.45:
        d = 2
    elif r < 0.65:
        d = 3
    else:
        d = min(max(1, int(rng.expovariate(1/6))), K)
    indices = rng.sample(range(K), min(d, K))
    result = bytearray(blocks[indices[0]])
    for idx in indices[1:]:
        for i in range(len(result)):
            result[i] ^= blocks[idx][i]
    return {'seed': seed, 'indices': indices, 'data': bytes(result)}

def fountain_encode(payload_bytes: bytes, block_size: int, num_droplets: int, seed0: int=12345):
    blocks = split_into_blocks(payload_bytes, block_size)
    droplets = []
    seed = seed0
    for i in range(num_droplets):
        droplet = fountain_generate_droplet(blocks, seed)
        droplets.append(droplet)
        seed += 1
    return blocks, droplets

def fountain_peeling_decode(droplets, block_size: int):
    D = [bytearray(d['data']) for d in droplets]
    idxs = [list(d['indices'])[:] for d in droplets]
    droplet_alive = [True]*len(droplets)
    maxidx = 0
    for indices in idxs:
        if indices:
            maxidx = max(maxidx, max(indices))
    K_guess = maxidx+1
    block_to_droplets = defaultdict(deque)
    for didx, indices in enumerate(idxs):
        for j in indices:
            block_to_droplets[j].append(didx)
    recovered = dict()
    queue = deque([d for d,i in enumerate(idxs) if len(i)==1])
    while queue:
        did = queue.popleft()
        if not droplet_alive[did]:
            continue
        if len(idxs[did]) != 1:
            continue
        blk = idxs[did][0]
        data = bytes(D[did])
        if blk in recovered:
            droplet_alive[did] = False
            continue
        recovered[blk] = data
        for other in list(block_to_droplets.get(blk, [])):
            if other == did or not droplet_alive[other]:
                continue
            for i in range(len(D[other])):
                D[other][i] ^= data[i]
            if blk in idxs[other]:
                idxs[other].remove(blk)
            if len(idxs[other]) == 1:
                queue.append(other)
        droplet_alive[did] = False
    if len(recovered) < K_guess:
        return None
    blocks = [recovered[i] for i in range(K_guess)]
    return blocks

# -------------------------
# Build oligos
# -------------------------
def build_oligos_from_payload_blocks(blocks, rs_codec, replicate=1, include_seed_in_index=True):
    parity_bytes = rs_codec.nsym
    oligos = []
    for idx, block in enumerate(blocks):
        bbits = ''.join(f"{b:08b}" for b in block)
        payload_dna = bits_to_dna_lookup(bbits)
        idx_bits = f"{idx:016b}"
        idx_dna = bits_to_dna_lookup(idx_bits)
        encoded = rs_codec.encode(block)
        ecc_bytes = encoded[-parity_bytes:]
        ecc_bits = ''.join(f"{b:08b}" for b in ecc_bytes)
        ecc_dna = bits_to_dna_lookup(ecc_bits)
        full = PRIMER_FWD + idx_dna + payload_dna + ecc_dna + PRIMER_REV
        for r in range(replicate):
            oligos.append((idx, r, full, {'type': 'block', 'orig_len': len(block)}))
    return oligos

def build_oligos_from_droplets(droplets, rs_codec, replicate=1):
    parity_bytes = rs_codec.nsym
    oligos = []
    for d in droplets:
        seed = d['seed']
        data = d['data']
        seed_hash = hashlib.sha256(str(seed).encode()).digest()
        idx_val = int.from_bytes(seed_hash[:2], 'big')
        idx_bits = f"{idx_val:016b}"
        idx_dna = bits_to_dna_lookup(idx_bits)
        bbits = ''.join(f"{b:08b}" for b in data)
        payload_dna = bits_to_dna_lookup(bbits)
        encoded = rs_codec.encode(data)
        ecc_bytes = encoded[-parity_bytes:]
        ecc_bits = ''.join(f"{b:08b}" for b in ecc_bytes)
        ecc_dna = bits_to_dna_lookup(ecc_bits)
        full = PRIMER_FWD + idx_dna + payload_dna + ecc_dna + PRIMER_REV
        for r in range(replicate):
            oligos.append((seed, r, full, {'type':'droplet','indices': d['indices']}))
    return oligos

# -------------------------
# Error simulation
# -------------------------
def simulate_errors_on_pool(oligo_pool, mutation_rate=0.0, dropout_rate=0.0, indel_rate=0.0):
    bases = ['A','T','C','G']
    survived = []
    for oid, rep, seq, info in oligo_pool:
        if random.random() < dropout_rate:
            continue
        seq_list = list(seq)
        i = 0
        out = []
        while i < len(seq_list):
            base = seq_list[i]
            if random.random() < indel_rate:
                if random.random() < 0.5:
                    i += 1
                    continue
                else:
                    out.append(random.choice([b for b in bases if b != base]))
            if random.random() < mutation_rate:
                out.append(random.choice([b for b in bases if b != base]))
            else:
                out.append(base)
            i += 1
        mutated = ''.join(out)
        survived.append((oid, rep, mutated, info))
    return survived

# -------------------------
# Decoding
# -------------------------
def strip_and_decode_inner_rs(seq, rs_codec):
    core = seq
    if core.startswith(PRIMER_FWD):
        core = core[len(PRIMER_FWD):]
    if core.endswith(PRIMER_REV):
        core = core[:-len(PRIMER_REV)]
    idx_dna = core[:16]
    payload_plus_ecc = core[16:]
    parity_b = rs_codec.nsym
    ecc_dna_len = (parity_b * 8) // 2
    if len(payload_plus_ecc) < ecc_dna_len:
        raise ValueError("oligo too short")
    payload_dna = payload_plus_ecc[:-ecc_dna_len]
    ecc_dna = payload_plus_ecc[-ecc_dna_len:]
    idx_bits = dna_to_bits_lookup(idx_dna)
    idx_bits = idx_bits[:16]
    idx_int = int(idx_bits, 2)
    payload_bits = dna_to_bits_lookup(payload_dna)
    if len(payload_bits) % 8 != 0:
        payload_bits = payload_bits[:-(len(payload_bits) % 8)]
    payload_bytes = bits_to_bytes_exact(payload_bits)
    ecc_bits = dna_to_bits_lookup(ecc_dna)
    if len(ecc_bits) % 8 != 0:
        ecc_bits = ecc_bits[:-(len(ecc_bits) % 8)]
    ecc_bytes = bits_to_bytes_exact(ecc_bits)
    decoded = rs_codec.decode(payload_bytes + ecc_bytes)[0]
    return idx_int, decoded

def try_decode_pool_rs_simple(corrupted_pool, rs_codec, mode='block'):
    success = 0
    fail = 0
    failed_ids = []
    recovered = {}
    droplets = []
    for oid, rep, seq, info in corrupted_pool:
        try:
            idx_int, decoded_bytes = strip_and_decode_inner_rs(seq, rs_codec)
            success += 1
            if info.get('type') == 'block' and mode=='block':
                recovered[idx_int] = decoded_bytes
            elif info.get('type') == 'droplet' and mode=='droplet':
                droplets.append({'seed': oid, 'indices': info.get('indices',[]), 'data': decoded_bytes})
            else:
                recovered[idx_int] = decoded_bytes
        except Exception:
            fail += 1
            failed_ids.append((oid, rep))
    if mode == 'block':
        return recovered, {'success': success, 'fail': fail, 'failed': failed_ids}
    else:
        return droplets, {'success': success, 'fail': fail, 'failed': failed_ids}

# -------------------------
# UI & orchestrator
# -------------------------
st.header("Simulator overview & controls")

st.write("This demo includes: Bytes→Bits→DNA, per-oligo RS inner-code, optional DNA-Fountain outer-code, replication, substitution+indel errors, RS decode, fountain peeling decode, success metrics and sweep plot.")

# Upload
uploaded = st.file_uploader("Choose a file", type=["txt","png","jpg","jpeg","mp3","wav","bin"])
if not uploaded:
    st.info("Upload a file to start. Keep small for quick testing.")
    st.stop()

raw = uploaded.read()
st.success(f"Uploaded: {uploaded.name} ({len(raw)} bytes)")

# Step A: Bytes -> Bits -> DNA
st.header("1) Bytes → Bits → DNA (ATCG)")
bit_string = bytes_to_bits(raw)
st.write(f"Total bytes: {len(raw)}, total bits: {len(bit_string)}")
st.code(bit_string[:400] + ("..." if len(bit_string) > 400 else ""))

dna_seq = bits_to_dna_lookup(bit_string)
st.write(f"Total DNA bases: {len(dna_seq)}")
st.code(dna_seq[:400] + ("..." if len(dna_seq) > 400 else ""))

# DNA stats
st.header("2) DNA statistics")
counts = base_counts(dna_seq)
df_counts = pd.DataFrame.from_dict(counts, orient='index', columns=['count'])
df_counts['fraction'] = (df_counts['count'] / len(dna_seq)).round(4)
st.table(df_counts)
st.write(f"GC content: **{gc_content(dna_seq):.2f}%**")
homopolys = find_homopolymers(dna_seq, min_run=6)
if homopolys:
    st.warning(f"Found {len(homopolys)} homopolymer runs (≥6). Showing up to 10.")
    st.table(pd.DataFrame(homopolys[:10], columns=['base','run_len','start','end']))
else:
    st.info("No long homopolymers detected (≥6).")

# Sidebar controls
st.sidebar.header(f"{COMPANY_NAME} — Encoder / Simulation settings")
nsym = st.sidebar.slider("RS parity bytes (nsym)", min_value=4, max_value=40, value=10, step=2)
rs_codec = make_rs(nsym)
oligo_len = st.sidebar.slider("Oligo payload length (bases)", min_value=60, max_value=300, value=120, step=10)
replicate = st.sidebar.slider("Replication factor (copies per oligo)", 1, 5, 1, 1)

use_fountain = st.sidebar.checkbox("Enable DNA-Fountain outer code (droplets)", value=False)
if use_fountain:
    block_size = st.sidebar.number_input("Fountain block size (bytes)", min_value=16, max_value=4096, value=256, step=16)
    droplets_count = st.sidebar.number_input("Number of fountain droplets to generate", min_value=1, max_value=10000, value=max(1, len(raw)//block_size + 20), step=1)

# Build oligo pool
st.header("3) Build oligo pool (inner RS per oligo + primers), optionally with Fountain")
if use_fountain:
    st.write("Using Fountain outer-code: splitting payload into blocks and generating droplets.")
    blocks = split_into_blocks(raw, block_size)
    st.write(f"Payload split into {len(blocks)} blocks of {block_size} bytes (last padded).")
    droplets = []
    seed0 = random.randint(1, 1<<30)
    for i in range(droplets_count):
        droplets.append(fountain_generate_droplet(blocks, seed0 + i))
    oligo_blocks = build_oligos_from_payload_blocks(blocks, rs_codec, replicate=replicate)
    oligo_droplets = build_oligos_from_droplets(droplets, rs_codec, replicate=replicate)
    oligo_pool = oligo_droplets
    st.write(f"Generated {len(droplets)} droplets → {len(oligo_droplets)} oligos (with replication).")
else:
    oligo_payloads = []
    for i in range(0, len(dna_seq), oligo_len):
        payload = dna_seq[i:i+oligo_len]
        oligo_payloads.append(payload)
    payload_blocks_bytes = []
    for p in oligo_payloads:
        pbits = dna_to_bits_lookup(p)
        if len(pbits) % 8 != 0:
            pbits = pbits.ljust((len(pbits)//8+1)*8, '0')
        pbytes = bits_to_bytes_exact(pbits)
        payload_blocks_bytes.append(pbytes)
    oligo_pool = build_oligos_from_payload_blocks(payload_blocks_bytes, rs_codec, replicate=replicate)
    st.write(f"Built {len(oligo_pool)} oligos (replicated = {replicate}).")

# preview
st.subheader("Oligo pool preview (first 12)")
preview_rows = []
for entry in oligo_pool[:12]:
    oid, rep, seq, info = entry
    core = seq
    if core.startswith(PRIMER_FWD):
        core = core[len(PRIMER_FWD):]
    if core.endswith(PRIMER_REV):
        core = core[:-len(PRIMER_REV)]
    payload_preview = core[16:16+min(60, len(core)-16)]
    preview_rows.append({'id': oid, 'rep': rep, 'type': info.get('type'), 'payload_preview': payload_preview})
st.table(pd.DataFrame(preview_rows))

st.download_button("Download pool as FASTA (oligos)", "\n".join([f">oligo_{o[0]}_rep{o[1]}\n{o[2]}" for o in oligo_pool]), file_name="oligo_pool.fasta")

# Simulation controls
st.header("4) Error simulation controls")
mutation_rate = st.slider("Substitution mutation rate (%)", 0, 20, 2, step=1)/100.0
indel_rate = st.slider("Indel (insertion/deletion) rate per base (%)", 0, 10, 0, step=1)/100.0
dropout_rate = st.slider("Dropout rate per oligo (%)", 0, 50, 0, step=1)/100.0

run_sim = st.button("Run simulation (apply errors & decode)")

# Sweep controls
st.header("Optional: Error-sweep experiment (mutation % sweep)")
do_sweep = st.checkbox("Run mutation sweep (vary substitution rate and show curve)", value=False)
sweep_min = st.number_input("Sweep min mutation %", min_value=0.0, max_value=10.0, value=0.0, step=0.5)
sweep_max = st.number_input("Sweep max mutation %", min_value=0.5, max_value=20.0, value=6.0, step=0.5)
sweep_steps = st.number_input("Sweep steps", min_value=3, max_value=25, value=7, step=1)
sweep_trials = st.number_input("Trials per sweep point (avg)", min_value=1, max_value=20, value=3, step=1)

# Run sim
if run_sim:
    corrupted = simulate_errors_on_pool(oligo_pool, mutation_rate=mutation_rate, dropout_rate=dropout_rate, indel_rate=indel_rate)
    st.write(f"Survived after dropout: {len(corrupted)} / {len(oligo_pool)}")
    if use_fountain:
        droplets_decoded, metrics = try_decode_pool_rs_simple(corrupted, rs_codec, mode='droplet')
        st.write("Inner RS decode metrics:", metrics)
        if droplets_decoded is None or len(droplets_decoded)==0:
            st.error("No droplets successfully decoded; fountain cannot run.")
            recovered_bytes = None
        else:
            blocks_recovered = fountain_peeling_decode(droplets_decoded, block_size)
            if blocks_recovered is None:
                st.error("Fountain peeling decoder failed to recover all blocks.")
                recovered_bytes = None
            else:
                recovered_bytes = b''.join(blocks_recovered)
                recovered_bytes = recovered_bytes[:len(raw)]
                st.success("Fountain + RS recovered the payload (subject to padding trimming).")
    else:
        recovered_map, metrics = try_decode_pool_rs_simple(corrupted, rs_codec, mode='block')
        st.write("Inner RS decode metrics:", metrics)
        if not recovered_map:
            st.error("No blocks decoded successfully.")
            recovered_bytes = None
        else:
            recovered_indices = sorted(recovered_map.keys())
            recovered_bytes = b''.join(recovered_map[i] for i in recovered_indices)
            recovered_bytes = recovered_bytes[:len(raw)]
            if len(recovered_indices) * oligo_len * 0.5 < len(raw):
                st.warning("Partial reconstruction: some indices missing—outer code or replication would be needed for full recovery.")
    if 'metrics' in locals():
        success = metrics.get('success',0)
        fail = metrics.get('fail',0)
        st.subheader("Recovery metrics (inner RS stage)")
        st.write(f"✅ inner RS successes: {success}")
        st.write(f"❌ inner RS fails: {fail}")
        total = success + fail
        success_rate = (success/ total*100) if total>0 else 0.0
        st.metric("Inner RS success rate", f"{success_rate:.2f}%")
        st.bar_chart(pd.DataFrame({"Status":["Success","Failed"], "Count":[success, fail]}).set_index("Status"))
        if success_rate==100:
            st.success("Perfect inner-RS decode of received oligos.")
        elif success_rate>=70:
            st.warning("Partial inner-RS decode.")
        else:
            st.error("Inner-RS decode poor — try increasing RS parity or lowering errors.")
    if recovered_bytes:
        st.subheader("Final recovered file")
        st.download_button("Download recovered file", recovered_bytes, file_name="recovered_"+uploaded.name)
    else:
        st.info("No full recovered file available with current settings. Try lowering error rates, increasing RS parity, using Fountain, or increasing replication.")

# Sweep experiment
if do_sweep:
    st.header("Error-sweep: mutation % vs recovery success")
    sweep_rates = [sweep_min + i*(sweep_max - sweep_min)/(sweep_steps-1) for i in range(sweep_steps)]
    results = []
    for rate in sweep_rates:
        avg_success_rate = 0.0
        for t in range(sweep_trials):
            corrupted = simulate_errors_on_pool(oligo_pool, mutation_rate=rate/100.0, dropout_rate=dropout_rate, indel_rate=indel_rate)
            _, metrics = try_decode_pool_rs_simple(corrupted, rs_codec, mode = 'droplet' if use_fountain else 'block')
            succ = metrics.get('success',0)
            total = metrics.get('success',0) + metrics.get('fail',0)
            sr = (succ/total*100) if total>0 else 0.0
            avg_success_rate += sr
        avg_success_rate /= sweep_trials
        results.append((rate, avg_success_rate))
    df_sweep = pd.DataFrame(results, columns=['mutation_%','inner_RS_success_%'])
    st.line_chart(df_sweep.set_index('mutation_%'))
    st.dataframe(df_sweep)

# Footer with branding + deploy notes
st.write("---")
left, right = st.columns([3,1])
with left:
    st.markdown(f"Built by **{COMPANY_NAME}** — DNA Storage Lab")
    st.markdown("For production-grade pipelines, integrate this simulator with laboratory workflows (oligo ordering, sequencing FASTQ input, and primer design).")
    st.markdown(f"[Visit {COMPANY_NAME}](https://codecellai.vercel.app/)")
with right:
    if os.path.exists(LOGO_PATH):
        st.image(LOGO_PATH, width=140)
    else:
        st.write("No logo found.")
st.info("Notes: This is a simulator for learning/prototyping. Fountain implementation here is a simplified LT-style toy for pedagogy. For production, follow Erlich & Zielinski (Science 2017) and use optimized robust soliton parameters.")
