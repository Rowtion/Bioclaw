#!/usr/bin/env python3
"""
TorchDrug: Drug Discovery ML Examples
======================================
PyTorch-based graph neural networks for molecules and proteins.
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Check for torchdrug availability
try:
    import torch
    import torchdrug
    from torchdrug import datasets, models, tasks, data
    TORCHDRUG_AVAILABLE = True
except ImportError:
    TORCHDRUG_AVAILABLE = False
    print("Warning: torchdrug not installed. Some examples will be skipped.")


# ==============================================================================
# Example 1: Molecular Property Prediction with GNN
# ==============================================================================

def example_molecular_property():
    """Predict molecular properties using graph neural networks."""
    print("=" * 60)
    print("Example 1: Molecular Property Prediction")
    print("=" * 60)
    
    if not TORCHDRUG_AVAILABLE:
        print("TorchDrug not available. Skipping example.")
        return None
    
    try:
        # Load BBBP dataset (Blood-Brain Barrier Penetration)
        dataset = datasets.BBBP(
            path="~/.torchdrug/BBBP/",
            verbose=1
        )
        
        print(f"Dataset loaded: {len(dataset)} molecules")
        print(f"Tasks: {dataset.tasks}")
        print(f"Node feature dim: {dataset.node_feature_dim}")
        print(f"Edge feature dim: {dataset.edge_feature_dim}")
        
        # Split dataset
        train_set, valid_set, test_set = dataset.split()
        print(f"Train: {len(train_set)}, Valid: {len(valid_set)}, Test: {len(test_set)}")
        
        # Define GNN model (GIN - Graph Isomorphism Network)
        model = models.GIN(
            input_dim=dataset.node_feature_dim,
            hidden_dims=[256, 256, 256],
            edge_input_dim=dataset.edge_feature_dim,
            batch_norm=True,
            readout="mean"
        )
        
        print(f"\nModel: GIN")
        print(f"  Hidden dims: {model.hidden_dims}")
        print(f"  Readout: mean")
        
        # Create property prediction task
        task = tasks.PropertyPrediction(
            model,
            task=dataset.tasks,
            criterion="bce",
            metric=["auroc", "auprc"],
            normalization=False
        )
        
        # Sample forward pass
        sample = train_set[0]
        print(f"\nSample molecule graph:")
        print(f"  Nodes: {sample['graph'].num_node}")
        print(f"  Edges: {sample['graph'].num_edge}")
        print(f"  Label: {sample['targets']}")
        
        return task, dataset
        
    except Exception as e:
        print(f"Error loading molecular data: {e}")
        return None


# ==============================================================================
# Example 2: Protein Structure Modeling
# ==============================================================================

def example_protein_modeling():
    """Work with protein structures using GNNs."""
    print("\n" + "=" * 60)
    print("Example 2: Protein Structure Modeling")
    print("=" * 60)
    
    if not TORCHDRUG_AVAILABLE:
        print("TorchDrug not available. Skipping example.")
        return None
    
    try:
        # Load protein dataset
        dataset = datasets.EnzymeCommission(
            path="~/.torchdrug/EnzymeCommission/",
            verbose=1
        )
        
        print(f"Dataset loaded: {len(dataset)} proteins")
        print(f"Tasks: {len(dataset.tasks)} EC classes")
        
        # Sample protein
        sample = dataset[0]
        protein = sample['graph']
        
        print(f"\nSample protein:")
        print(f"  Residues (nodes): {protein.num_node}")
        print(f"  Edges: {protein.num_edge}")
        
        # Demonstrate GearNet model for proteins
        model = models.GearNet(
            input_dim=dataset.node_feature_dim,
            hidden_dims=[512, 512, 512],
            num_relation=8,
            edge_input_dim=dataset.edge_feature_dim,
            batch_norm=True,
            readout="mean"
        )
        
        print(f"\nModel: GearNet (protein-specific GNN)")
        print(f"  Hidden dims: {model.hidden_dims}")
        
        return model, dataset
        
    except Exception as e:
        print(f"Error loading protein data: {e}")
        return None


# ==============================================================================
# Example 3: Working with Molecules and RDKit
# ==============================================================================

def example_molecule_rdkit():
    """Convert between TorchDrug and RDKit molecules."""
    print("\n" + "=" * 60)
    print("Example 3: Molecule-RDKit Conversion")
    print("=" * 60)
    
    if not TORCHDRUG_AVAILABLE:
        print("TorchDrug not available. Skipping example.")
        return None
    
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, Descriptors
        
        # Create molecule from SMILES
        smiles_list = [
            "CCO",  # Ethanol
            "CC(=O)O",  # Acetic acid
            "c1ccccc1",  # Benzene
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
        ]
        
        print("Processing molecules from SMILES:")
        molecules = []
        
        for smiles in smiles_list:
            # SMILES to TorchDrug molecule
            mol = data.Molecule.from_smiles(smiles, atom_feature="default")
            molecules.append(mol)
            
            # Get properties
            rdkit_mol = Chem.MolFromSmiles(smiles)
            mw = Descriptors.MolecularWeight(rdkit_mol)
            logp = Descriptors.MolLogP(rdkit_mol)
            
            print(f"\n  SMILES: {smiles}")
            print(f"  Atoms: {mol.num_node}, Bonds: {mol.num_edge}")
            print(f"  MW: {mw:.2f}, LogP: {logp:.2f}")
        
        # Draw molecules
        rdkit_mols = [Chem.MolFromSmiles(s) for s in smiles_list]
        img = Draw.MolsToGridImage(rdkit_mols, molsPerRow=2, subImgSize=(300, 300))
        
        # Save if possible
        try:
            img.save('torchdrug_molecules.png')
            print("\nSaved: torchdrug_molecules.png")
        except:
            print("\nCould not save molecule image")
        
        return molecules
        
    except ImportError:
        print("RDKit not available. Skipping RDKit integration.")
        return None
    except Exception as e:
        print(f"Error: {e}")
        return None


# ==============================================================================
# Example 4: Knowledge Graph Embedding
# ==============================================================================

def example_knowledge_graph():
    """Knowledge graph completion for drug repurposing."""
    print("\n" + "=" * 60)
    print("Example 4: Knowledge Graph Completion")
    print("=" * 60)
    
    if not TORCHDRUG_AVAILABLE:
        print("TorchDrug not available. Skipping example.")
        return None
    
    try:
        # Try to load Hetionet (biomedical knowledge graph)
        dataset = datasets.Hetionet(
            path="~/.torchdrug/Hetionet/",
            verbose=1
        )
        
        print(f"Knowledge graph loaded:")
        print(f"  Entities: {dataset.num_entity}")
        print(f"  Relations: {dataset.num_relation}")
        print(f"  Triples: {len(dataset)}")
        
        # Define RotatE model for KG embedding
        model = models.RotatE(
            num_entity=dataset.num_entity,
            num_relation=dataset.num_relation,
            embedding_dim=128,
            gamma=12
        )
        
        # Create knowledge graph completion task
        task = tasks.KnowledgeGraphCompletion(
            model,
            criterion="bce",
            metric=["mr", "mrr", "hits@1", "hits@3", "hits@10"]
        )
        
        print(f"\nModel: RotatE")
        print(f"  Embedding dim: 128")
        print(f"  Evaluation metrics: MRR, Hits@K")
        
        return task, dataset
        
    except Exception as e:
        print(f"Error loading knowledge graph: {e}")
        print("Note: Hetionet is a large dataset that may need to be downloaded.")
        return None


# ==============================================================================
# Example 5: Retrosynthesis Planning
# ==============================================================================

def example_retrosynthesis():
    """Chemical retrosynthesis planning."""
    print("\n" + "=" * 60)
    print("Example 5: Retrosynthesis Planning")
    print("=" * 60)
    
    if not TORCHDRUG_AVAILABLE:
        print("TorchDrug not available. Skipping example.")
        return None
    
    try:
        # Load USPTO reaction dataset
        dataset = datasets.USPTO50k(
            path="~/.torchdrug/USPTO50k/",
            verbose=1
        )
        
        print(f"Reaction dataset loaded: {len(dataset)} reactions")
        
        # Sample reaction
        sample = dataset[0]
        print(f"\nSample reaction:")
        print(f"  Type: {sample['reaction']}")
        print(f"  Center identified: {sample.get('center', 'N/A')}")
        
        # Retrosynthesis models
        print(f"\nRetrosynthesis pipeline:")
        print(f"  1. CenterIdentification: Predict reaction center")
        print(f"  2. SynthonCompletion: Complete reactants")
        print(f"  3. Combine into end-to-end pipeline")
        
        # Define center identification model
        center_model = models.RGCN(
            input_dim=dataset.node_feature_dim,
            hidden_dims=[256, 256],
            num_relation=dataset.num_relation,
            batch_norm=True
        )
        
        center_task = tasks.CenterIdentification(
            center_model,
            feature=dataset.node_feature_dim
        )
        
        print(f"\nCenter identification model: RGCN")
        
        return center_task, dataset
        
    except Exception as e:
        print(f"Error loading retrosynthesis data: {e}")
        return None


# ==============================================================================
# Example 6: Molecular Generation
# ==============================================================================

def example_molecular_generation():
    """Generate novel molecular structures."""
    print("\n" + "=" * 60)
    print("Example 6: Molecular Generation")
    print("=" * 60)
    
    if not TORCHDRUG_AVAILABLE:
        print("TorchDrug not available. Skipping example.")
        return None
    
    try:
        print("Molecular Generation Approaches:")
        print("  1. GCPN (Graph Convolutional Policy Network)")
        print("     - RL-based molecular optimization")
        print("     - Property-guided generation")
        
        print("\n  2. GraphAutoregressiveFlow")
        print("     - Flow-based generation")
        print("     - Normalizing flows for molecules")
        
        print("\n  3. Autoregressive Generation")
        print("     - Sequential atom/bond addition")
        print("     - SMILES-based generation")
        
        # Define GCPN model for demonstration
        print("\nExample GCPN configuration:")
        print("  - Node features: atom type, charge, etc.")
        print("  - Edge features: bond type")
        print("  - Policy network: GNN + MLP")
        print("  - Reward: molecular properties (QED, LogP, etc.)")
        
        return None
        
    except Exception as e:
        print(f"Error: {e}")
        return None


# ==============================================================================
# Example 7: Working with Graph Data
# ==============================================================================

def example_graph_data():
    """Understanding TorchDrug graph data structures."""
    print("\n" + "=" * 60)
    print("Example 7: Graph Data Structures")
    print("=" * 60)
    
    if not TORCHDRUG_AVAILABLE:
        print("TorchDrug not available. Skipping example.")
        return None
    
    try:
        # Create a simple molecule
        smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"  # Ibuprofen
        mol = data.Molecule.from_smiles(smiles)
        
        print(f"Molecule: {smiles}")
        print(f"  Num atoms: {mol.num_node}")
        print(f"  Num bonds: {mol.num_edge}")
        print(f"  Atom types: {mol.atom_type}")
        print(f"  Bond types: {mol.edge_list[:, 2]}")
        
        # Graph attributes
        print(f"\nGraph attributes:")
        print(f"  Node features shape: {mol.node_feature.shape}")
        print(f"  Edge features shape: {mol.edge_feature.shape}" if hasattr(mol, 'edge_feature') else "  No edge features")
        
        # Convert to PyG format (if needed)
        print(f"\nGraph conversion:")
        print(f"  Can convert to PyTorch Geometric Data object")
        print(f"  Can convert to NetworkX graph")
        print(f"  Can convert to DGL graph")
        
        # Batch processing
        print(f"\nBatch processing:")
        print(f"  Use data.PackedGraph for variable-size graphs")
        print(f"  Automatic padding and masking")
        
        return mol
        
    except Exception as e:
        print(f"Error: {e}")
        return None


# ==============================================================================
# Example 8: Custom Dataset and Data Loading
# ==============================================================================

def example_custom_dataset():
    """Create custom datasets for molecular data."""
    print("\n" + "=" * 60)
    print("Example 8: Custom Dataset Creation")
    print("=" * 60)
    
    if not TORCHDRUG_AVAILABLE:
        print("TorchDrug not available. Skipping example.")
        return None
    
    print("Creating custom dataset:")
    print("""
from torchdrug import data

class CustomDataset(data.Dataset):
    def __init__(self, smiles_list, labels):
        super().__init__()
        self.smiles_list = smiles_list
        self.labels = labels
        
    def get_item(self, index):
        smiles = self.smiles_list[index]
        label = self.labels[index]
        
        # Convert SMILES to graph
        graph = data.Molecule.from_smiles(smiles)
        
        return {
            "graph": graph,
            "label": label
        }
    
    def __len__(self):
        return len(self.smiles_list)
""")
    
    # Demonstrate with sample data
    smiles_data = [
        ("CCO", 1),  # Ethanol - class 1
        ("CC(=O)O", 1),  # Acetic acid - class 1
        ("CCCC", 0),  # Butane - class 0
        ("CC(C)C", 0),  # Isobutane - class 0
    ]
    
    print(f"\nSample data: {len(smiles_data)} molecules")
    for smiles, label in smiles_data:
        mol = data.Molecule.from_smiles(smiles)
        print(f"  {smiles}: {mol.num_node} atoms, label={label}")
    
    print(f"\nData loading:")
    print(f"  - Use torch.utils.data.DataLoader")
    print(f"  - Set collate_fn for graph batching")
    print(f"  - Use scaffold split for molecular datasets")
    
    return smiles_data


# ==============================================================================
# Main Execution
# ==============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("TorchDrug: Drug Discovery ML Examples")
    print("=" * 70)
    
    if not TORCHDRUG_AVAILABLE:
        print("\nNote: TorchDrug is not installed.")
        print("To install: pip install torchdrug")
        print("Examples will run in demonstration mode.\n")
    
    try:
        example_molecular_property()
        example_protein_modeling()
        example_molecule_rdkit()
        example_knowledge_graph()
        example_retrosynthesis()
        example_molecular_generation()
        example_graph_data()
        example_custom_dataset()
        
        print("\n" + "=" * 70)
        print("All examples completed!")
        print("=" * 70)
        
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
