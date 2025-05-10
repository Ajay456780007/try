import os
import gzip
from Bio import SeqIO
import pandas as pd
from urllib.parse import unquote

# === Paths ===
dna_dir = "dataset/dataset1/dna_chromosomes"
gff3_dir = "dataset/dataset1/gff3_files"

# === Collect all FASTA and GFF3 files ===
fasta_files = sorted([
    os.path.join(dna_dir, f) for f in os.listdir(dna_dir)
    if f.lower().endswith(".fa.gz")
])
gff3_files = sorted([
    os.path.join(gff3_dir, f) for f in os.listdir(gff3_dir)
    if f.lower().endswith(".gff3")
])

# === Parse GFF3 attributes ===
def parse_attributes(attr_str):
    attr_dict = {}
    for pair in attr_str.strip().split(";"):
        if "=" in pair:
            key, value = pair.split("=", 1)
            attr_dict[key.strip()] = unquote(value.strip())
    return attr_dict

# === Function to parse GFF3 and extract gene entries for a given chromosome ===
def parse_gff3(gff3_file, chrom_id):
    genes = []
    with open(gff3_file, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "gene":
                continue
            if parts[0] != chrom_id:
                continue  # skip if this line refers to a different chromosome

            start = int(parts[3]) - 1  # Convert to 0-based index
            end = int(parts[4])
            strand = parts[6]
            attrs = parse_attributes(parts[8])
            gene_id = attrs.get("ID", "NA")

            genes.append((start, end, strand, gene_id))
    return genes

# === Extract gene sequences ===
gene_sequences = []

for fasta_path, gff_path in zip(fasta_files, gff3_files):
    print(f"Processing: {os.path.basename(fasta_path)} with {os.path.basename(gff_path)}")

    # Read the chromosome sequence (support gzipped or uncompressed)
    if fasta_path.endswith(".gz"):
        with gzip.open(fasta_path, "rt", encoding="utf-8") as f:
            record = next(SeqIO.parse(f, "fasta"))
    else:
        with open(fasta_path, "r", encoding="utf-8") as f:
            record = next(SeqIO.parse(f, "fasta"))

    chrom_seq = record.seq
    chrom_id = record.id

    # Parse gene entries from the GFF3 file
    genes = parse_gff3(gff_path, chrom_id)

    for start, end, strand, gene_id in genes:
        # Boundary check to avoid errors
        if start < 0 or end > len(chrom_seq):
            continue

        gene_seq = chrom_seq[start:end]
        if strand == "-":
            gene_seq = gene_seq.reverse_complement()

        gene_sequences.append({
            "gene_id": gene_id,
            "chrom": chrom_id,
            "start": start,
            "end": end,
            "strand": strand,
            "sequence": str(gene_seq)
        })

# === Convert to DataFrame ===
df_genes = pd.DataFrame(gene_sequences)

# === Preview first few entries ===
print(df_genes.head())
geo_path="dataset/dataset1/geo_files/genes_to_alias_ids.tsv"
df = pd.read_csv(geo_path, sep='\t')
alias_path="dataset/dataset1/expression level TPM/abundance.tsv"
df_alias = pd.read_csv(alias_path, sep='\t')
df.rename(columns={"Zm00001eb000010":"id1","B73 Zm00001eb.1":"id2","Zm00001d027230":"gene_alias_id"})

import pandas as pd

# Fix column names in 'df'
df.columns = ['id1', 'id2', 'gene_alias_id', 'AGPv4_Zm00001d.2']

# Step 1: Clean the 'gene_id' column in df_genes
df_genes['gene_id_clean'] = df_genes['gene_id'].str.replace('gene:', '', regex=False)

# Step 2: Create a mapping from 'id1' to 'gene_alias_id'
id_to_alias = df.set_index('id1')['gene_alias_id'].to_dict()

# Step 3: Map gene_alias_id to df_genes
df_genes['alias_id'] = df_genes['gene_id_clean'].map(id_to_alias)

# Step 4: Clean up the temporary column
df_genes.drop(columns=['gene_id_clean'], inplace=True)

# Check the results
print(df_genes[['gene_id', 'alias_id']].head(50))


# Drop rows where alias_id is NaN
df_genes = df_genes.dropna(subset=['alias_id'])

# Reset index if needed
df_genes = df_genes.reset_index(drop=True)

# Optional: check that it's gone
print(df_genes['alias_id'].isna().sum())



import pandas as pd

# Step 1: Clean up df_alias to remove transcript suffix
df_alias['clean_id'] = df_alias['target_id'].str.replace(r'_T\d+$', '', regex=True)

# Step 2: Group by clean_id and sum or average TPM if needed (in case multiple transcripts per gene)
# Here we'll sum TPMs for all isoforms of a gene
tpm_by_gene = df_alias.groupby('clean_id')['tpm'].sum().reset_index()

# Step 3: Merge df_genes with this TPM info using alias_id == clean_id
df_genes = df_genes.merge(tpm_by_gene, how='left', left_on='alias_id', right_on='clean_id')

# Step 4: Rename the column to tpm_value and drop clean_id
df_genes = df_genes.rename(columns={'tpm': 'tpm_value'}).drop(columns=['clean_id'])

# Done
print(df_genes.head())


# Mapping dictionary
base_map = {'A': 0, 'T': 1, 'G': 2, 'C': 3}

# Function to encode a DNA sequence string
def encode_sequence(seq):
    return [base_map.get(base, -1) for base in seq.upper()]  # -1 for unknown bases like N

# Apply to each row in df_genes['sequence']
df_genes['encoded_sequence'] = df_genes['sequence'].apply(encode_sequence)

# Done
print(df_genes[['sequence', 'encoded_sequence']].head())


df_genes.drop("sequence",axis=1,inplace=True)


max_len = df_genes['encoded_sequence'].apply(len).max()
print("Maximum encoded sequence length:", max_len)




# Step 1: Store lengths of all encoded sequences
sequence_lengths = df_genes['encoded_sequence'].apply(len)

# Step 2: Count how many are greater than 50,000
num_greater_than_40000 = (sequence_lengths > 40000).sum()

# Output results
print("Total sequences:", len(sequence_lengths))
print("Sequences > 50,000 bases:", num_greater_than_40000)



# Fixed length
FIXED_LEN = 500
PAD_VALUE = 0  # A = 0

def pad_or_truncate(seq):
    if len(seq) > FIXED_LEN:
        return seq[:FIXED_LEN]  
    else:
        return seq + [PAD_VALUE] * (FIXED_LEN - len(seq))  # pad

# Apply to each sequence
df_genes['encoded_50k'] = df_genes['encoded_sequence'].apply(pad_or_truncate)

# Check shape of one example
print(len(df_genes['encoded_50k'].iloc[0])) 
  
  
print(df_genes.head())

df_genes['strand'] = df_genes['strand'].map({'+': 1, '-': 0})

df_genes.rename(columns={"encoded_50k":"sequence"},inplace=True)

df_genes["sequence"]

import numpy as np
def fix_sequence(seq):
    return [4 if val == -1 else val for val in seq]

df_genes['sequence'] = df_genes['sequence'].apply(fix_sequence)
x = np.array(df_genes['sequence'].tolist(), dtype=np.uint8)

from sklearn.model_selection import train_test_split

y=df_genes["tpm_value"]


print("this is x shape:",x.shape)
print("this is y shape:",y.shape) 

import numpy as np
y=np.array(y)


from sklearn.model_selection import train_test_split
x_train,x_test,y_train,y_test=train_test_split(x,y,test_size=0.25,random_state=25)

print(x_train.shape)
print(y_train.shape)
print(x_test.shape)
print(y_test.shape)

import os
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

import numpy as np
import tensorflow as tf
from tensorflow.keras.utils import Sequence
from tensorflow.keras import layers
import math

# Enable float16 precision for speed and memory savings
from tensorflow.keras import mixed_precision
mixed_precision.set_global_policy('mixed_float16')

# ------------ Step 1: Load and Trim the Dataset ------------

# Reduce sequence length and sample count for faster CPU training
x = np.load("x.npy", mmap_mode='r')[:2000, :2000]
y = np.load("y.npy", mmap_mode='r')[:2000]

np.save("x_trimmed.npy", x)
np.save("y_trimmed.npy", y)

# ------------ Step 2: Data Generator ------------
import numpy as np

class DNADataGenerator(tf.keras.utils.Sequence):
    def __init__(self, dna_data, labels, batch_size):
        self.dna_data = dna_data
        self.labels = labels
        self.batch_size = batch_size
        self.start = start
        self.end = len(self.x) if end is None else end
        self.indexes = np.arange(self.start, self.end)

    def __len__(self):
        return int(np.ceil((self.end - self.start) / self.batch_size))

    def __getitem__(self, index):
        start = self.start + index * self.batch_size
        end = min(start + self.batch_size, self.end)

        batch_x = self.x[start:end].astype(np.int32)
        batch_y = self.y[start:end].astype(np.float32)

        return batch_x, batch_y





# ------------ Step 3: Train/Test Split ------------
split_idx = int(2000 * 0.75)
# # Use the updated generator with only DNA and target data
train_gen = DNADataGenerator("x_trimmed.npy", "y_trimmed.npy", batch_size=4, start=0, end=split_idx)
test_gen = DNADataGenerator("x_trimmed.npy", "y_trimmed.npy", batch_size=4, start=split_idx)


# ------------ Step 4: Positional and Noise Embeddings ------------

def get_positional_encoding(seq_len, model_dim):
    angle_rads = np.arange(seq_len)[:, np.newaxis] / np.power(
        10000, (2 * (np.arange(model_dim)[np.newaxis, :] // 2)) / np.float32(model_dim)
    )
    angle_rads[:, 0::2] = np.sin(angle_rads[:, 0::2])
    angle_rads[:, 1::2] = np.cos(angle_rads[:, 1::2])
    return tf.constant(angle_rads[np.newaxis, ...], dtype=tf.float32)

class SinusoidalEmbedding(layers.Layer):
    def __init__(self, model_dim):
        super().__init__()
        self.model_dim = model_dim

    def call(self, x):
        half_dim = self.model_dim // 2
        freqs = tf.exp(tf.linspace(tf.math.log(1.0), tf.math.log(1000.0), half_dim))
        x = tf.cast(x, tf.float32)
        angles = 2.0 * math.pi * x * freqs
        sin = tf.sin(angles)
        cos = tf.cos(angles)
        return tf.concat([sin, cos], axis=-1)[..., tf.newaxis]

# ------------ Step 5: Transformer Block ------------

class TransformerBlock(layers.Layer):
    def __init__(self, model_dim, heads, ff_dim, rate=0.1):
        super().__init__()
        self.att = layers.MultiHeadAttention(num_heads=heads, key_dim=model_dim)
        self.ffn = tf.keras.Sequential([
            layers.Dense(ff_dim, activation='gelu'),
            layers.Dense(model_dim),
        ])
        self.norm1 = layers.LayerNormalization(epsilon=1e-6)
        self.norm2 = layers.LayerNormalization(epsilon=1e-6)
        self.dropout1 = layers.Dropout(rate)
        self.dropout2 = layers.Dropout(rate)

    def call(self, x, training=False):
        attn_output = self.att(x, x)
        out1 = self.norm1(x + self.dropout1(attn_output, training=training))
        ffn_output = self.ffn(out1)
        return self.norm2(out1 + self.dropout2(ffn_output, training=training))

# ------------ Step 6: Full Model Definition ------------

class GeneExpressionTransformer(tf.keras.Model):
    def __init__(self, seq_len=2000, model_dim=64, num_heads=2, ff_dim=128, depth=2):
        super().__init__()
        self.seq_len = seq_len
        self.model_dim = model_dim

        self.embedding = tf.keras.layers.Embedding(input_dim=5, output_dim=model_dim)
        self.pos_encoding = tf.cast(get_positional_encoding(seq_len, model_dim), dtype=tf.float16)
        self.time_emb = SinusoidalEmbedding(model_dim)

        self.transformer_blocks = [
            TransformerBlock(model_dim, num_heads, ff_dim) for _ in range(depth)
        ]

        self.global_pool = tf.keras.layers.GlobalAveragePooling1D()
        self.output_dense = tf.keras.layers.Dense(1, dtype='float32')  # Output back to float32

    def call(self, inputs, training=False):
        dna_seq = inputs
        noise_var = tf.zeros((tf.shape(dna_seq)[0], 1))

        x = self.embedding(dna_seq)
        x += self.pos_encoding[:, :tf.shape(x)[1], :]

        noise_emb = self.time_emb(noise_var)
        noise_emb = tf.transpose(noise_emb, [0, 2, 1])
        noise_emb = tf.tile(noise_emb, [1, tf.shape(x)[1], 1])
        noise_emb = tf.cast(noise_emb, dtype=tf.float16)  # 🔧 Cast to float16 to match x
        x += noise_emb

        for block in self.transformer_blocks:
            x = block(x, training=training)

        x = self.global_pool(x)
        return self.output_dense(x)

# ------------ Step 7: Compile and Train ------------

import tensorflow as tf
from tensorflow.keras import layers, models

# Define the Diffusion Transformer Model (assuming it's implemented already)
class DiffusionTransformer(tf.keras.Model):
    def __init__(self, seq_len=4000, model_dim=128, num_heads=4, ff_dim=256, depth=4):
        super(DiffusionTransformer, self).__init__()
        self.seq_len = seq_len
        self.model_dim = model_dim

        # Diffusion transformer architecture here...
        self.embedding = layers.Embedding(input_dim=5, output_dim=model_dim)
        self.pos_encoding = get_positional_encoding(seq_len, model_dim)
        self.time_emb = SinusoidalEmbedding(model_dim)
        self.transformer_blocks = [
            TransformerBlock(model_dim, num_heads, ff_dim) for _ in range(depth)
        ]
        self.global_pool = layers.GlobalAveragePooling1D()
        self.output_dense = layers.Dense(1)

    def call(self, inputs, training=False):
        dna_seq = inputs
        noise_var = tf.zeros((tf.shape(dna_seq)[0], 1))

        x = self.embedding(inputs)  # Use 'inputs', not 'x'
        x = tf.cast(x, dtype=tf.float32)
        x += self.pos_encoding[:, :tf.shape(x)[1], :]

        noise_emb = self.time_emb(noise_var)
        noise_emb = tf.transpose(noise_emb, [0, 2, 1])
        noise_emb = tf.tile(noise_emb, [1, tf.shape(x)[1], 1])
        x += noise_emb

        for block in self.transformer_blocks:
            x = block(x, training=training)

        x = self.global_pool(x)
        return self.output_dense(x)

# Define HGNN (Hybrid Graph Neural Network)
class HGNNConv(tf.keras.layers.Layer):
    def __init__(self, input_dim, output_dim):
        super(HGNNConv, self).__init__()
        self.weight = self.add_weight(shape=(input_dim, output_dim),
                                      initializer='glorot_uniform',
                                      trainable=True)
        self.bias = self.add_weight(shape=(output_dim,),
                                    initializer='zeros',
                                    trainable=True)

    def call(self, x, G):
        x = tf.matmul(x, self.weight) + self.bias
        x = tf.matmul(G, x)
        return x


class HGNNEmbedding(tf.keras.Model):
    def __init__(self, input_dim, hidden_dim, dropout_rate=0.5):
        super(HGNNEmbedding, self).__init__()
        self.hgc1 = HGNNConv(input_dim, hidden_dim)
        self.hgc2 = HGNNConv(hidden_dim, hidden_dim)
        self.dropout = dropout_rate

    def call(self, x, G, training=False):
        x = tf.nn.relu(self.hgc1(x, G))
        if training:
            x = tf.nn.dropout(x, rate=self.dropout)
        x = tf.nn.relu(self.hgc2(x, G))
        return x

# Define Hybrid Attention Fusion Model
class HybridAttentionFusion(tf.keras.Model):
    def __init__(self, vocab_size=5, embed_dim=64, conv_channels=64, kernel_sizes=[7, 5, 3], dropout_rate=0.2):
        super(HybridAttentionFusion, self).__init__()

        self.embedding = layers.Embedding(input_dim=vocab_size, output_dim=embed_dim, mask_zero=True)

        # Convolutional layers for DNA
        self.conv1 = layers.Conv1D(conv_channels, kernel_sizes[0], activation='relu')
        self.conv2 = layers.Conv1D(conv_channels * 2, kernel_sizes[1], activation='relu')
        self.conv3 = layers.Conv1D(conv_channels * 4, kernel_sizes[2], activation='relu')
        self.global_pool = layers.GlobalMaxPooling1D()

        # Fully connected layers
        self.dropout = layers.Dropout(dropout_rate)
        self.fc1 = layers.Dense(1024, activation='relu')
        self.fc2 = layers.Dense(512, activation='relu')
        self.out = layers.Dense(1, activation='sigmoid')  # or use softmax for multi-class

    def call(self, dna_input):
        # Embedding
        dna = self.embedding(dna_input)
        dna = self.conv1(dna)
        dna = self.conv2(dna)
        dna = self.conv3(dna)

        # Global pooling
        dna_feat = self.global_pool(dna)

        # Dense layers
        x = self.dropout(dna_feat)
        x = self.fc1(x)
        x = self.dropout(x)
        x = self.fc2(x)
        return self.out(x)


# Combined Model
class CombinedModel(tf.keras.Model):
    def __init__(self, seq_len, vocab_size, embed_dim, conv_channels, kernel_sizes, hidden_dim, num_classes):
        super(CombinedModel, self).__init__()

        # Instantiate individual models
        self.diffusion_transformer = DiffusionTransformer(seq_len=seq_len)
        self.hgnn_embedding = HGNNEmbedding(input_dim=5, hidden_dim=hidden_dim)
        self.hybrid_attention_fusion = HybridAttentionFusion(vocab_size, embed_dim, conv_channels, kernel_sizes)

        self.final_fc = layers.Dense(num_classes, activation='sigmoid')

    def call(self, dna_input, training=False):
        # Get diffusion transformer output
        diffusion_output = self.diffusion_transformer(dna_input, training=training)

        # Get hybrid attention fusion output (without chip_input)
        hybrid_attention_output = self.hybrid_attention_fusion(dna_input)

        # Concatenate all outputs
        combined_output = tf.concat([diffusion_output, hybrid_attention_output], axis=1)

        # Final fully connected layer
        return self.final_fc(combined_output)


# Model initialization
seq_len = 2000  # Length of the DNA sequence
vocab_size = 5  # A typical vocabulary size for DNA sequences (A, C, G, T, padding)
embed_dim = 64
conv_channels = 64
kernel_sizes = [7, 5, 3]
hidden_dim = 128
num_classes = 1  # For regression (use 1), for classification change accordingly

model = CombinedModel(seq_len, vocab_size, embed_dim, conv_channels, kernel_sizes, hidden_dim, num_classes)

# Compile and fit the model
model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=1e-4), loss='mse', metrics=['mae'])
model.fit(train_gen, validation_data=test_gen, epochs=5)

