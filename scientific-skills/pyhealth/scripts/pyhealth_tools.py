#!/usr/bin/env python3
"""
PyHealth Tools - Healthcare AI toolkit utilities
Clinical prediction models, medical coding, and dataset management
"""

import argparse
import json
import numpy as np
from typing import List, Dict, Optional


def list_available_datasets() -> List[Dict]:
    """List available healthcare datasets in PyHealth."""
    datasets = [
        {"name": "MIMIC3Dataset", "description": "MIMIC-III critical care database", "type": "EHR"},
        {"name": "MIMIC4Dataset", "description": "MIMIC-IV critical care database", "type": "EHR"},
        {"name": "eICUDataset", "description": "eICU Collaborative Research Database", "type": "EHR"},
        {"name": "OMOPDataset", "description": "OMOP Common Data Model", "type": "EHR"},
        {"name": "SleepEDFDataset", "description": "Sleep-EDF polysomnography", "type": "Signal"},
        {"name": "SHHSDataset", "description": "Sleep Heart Health Study", "type": "Signal"},
        {"name": "ISRUCDataset", "description": "ISRUC sleep dataset", "type": "Signal"},
        {"name": "COVID19CXRDataset", "description": "COVID-19 Chest X-ray", "type": "Image"},
    ]
    return datasets


def list_prediction_tasks() -> List[Dict]:
    """List available clinical prediction tasks."""
    tasks = [
        {"name": "mortality_prediction", "description": "In-hospital mortality prediction", "mode": "binary"},
        {"name": "readmission_prediction", "description": "30-day readmission prediction", "mode": "binary"},
        {"name": "length_of_stay_prediction", "description": "Length of stay prediction", "mode": "multiclass"},
        {"name": "drug_recommendation", "description": "Medication recommendation", "mode": "multilabel"},
        {"name": "sleep_staging", "description": "Sleep stage classification", "mode": "multiclass"},
        {"name": "medical_code_prediction", "description": "ICD code prediction from text", "mode": "multilabel"},
    ]
    return tasks


def list_models() -> List[Dict]:
    """List available models in PyHealth."""
    models = [
        {"name": "LogisticRegression", "type": "Baseline", "description": "Logistic regression classifier"},
        {"name": "MLP", "type": "Deep Learning", "description": "Multi-layer perceptron"},
        {"name": "CNN", "type": "Deep Learning", "description": "Convolutional neural network"},
        {"name": "RNN", "type": "Deep Learning", "description": "Recurrent neural network (LSTM/GRU)"},
        {"name": "Transformer", "type": "Deep Learning", "description": "Transformer architecture"},
        {"name": "RETAIN", "type": "Healthcare", "description": "REverse Time AttentIoN network"},
        {"name": "SafeDrug", "type": "Healthcare", "description": "Safe drug recommendation with DDI"},
        {"name": "GAMENet", "type": "Healthcare", "description": "Graph augmented MEmory network"},
        {"name": "StageNet", "type": "Healthcare", "description": "Stage-aware neural network"},
        {"name": "AdaCare", "type": "Healthcare", "description": "Adaptive healthcare model"},
    ]
    return models


def generate_mortality_prediction_template(output_path: str):
    """Generate a template script for mortality prediction."""
    template = '''#!/usr/bin/env python3
"""
Mortality Prediction using PyHealth
Example: MIMIC-IV mortality prediction with RETAIN model
"""

from pyhealth.datasets import MIMIC4Dataset
from pyhealth.tasks import mortality_prediction_mimic4_fn
from pyhealth.datasets import split_by_patient, get_dataloader
from pyhealth.models import RETAIN
from pyhealth.trainer import Trainer

# 1. Load dataset
print("Loading MIMIC-IV dataset...")
dataset = MIMIC4Dataset(
    root="/path/to/mimic4",
    tables=["diagnoses", "procedures", "prescriptions"],
    code_mapping={},
    dev=True  # Use dev mode for testing
)
print(dataset.stats())

# 2. Set task
print("Setting mortality prediction task...")
sample_dataset = dataset.set_task(mortality_prediction_mimic4_fn)
print(f"Samples: {len(sample_dataset)}")

# 3. Split data
print("Splitting data...")
train_dataset, val_dataset, test_dataset = split_by_patient(
    sample_dataset, ratios=[0.7, 0.1, 0.2], seed=42
)

# 4. Create dataloaders
train_loader = get_dataloader(train_dataset, batch_size=64, shuffle=True)
val_loader = get_dataloader(val_dataset, batch_size=64)
test_loader = get_dataloader(test_dataset, batch_size=64)

# 5. Initialize model
print("Initializing RETAIN model...")
model = RETAIN(
    dataset=sample_dataset,
    feature_keys=["conditions", "procedures", "drugs"],
    label_key="label",
    mode="binary",
    embedding_dim=128,
    hidden_dim=128
)

# 6. Train
print("Training model...")
trainer = Trainer(model=model, device="cuda")
trainer.train(
    train_dataloader=train_loader,
    val_dataloader=val_loader,
    epochs=50,
    optimizer="Adam",
    learning_rate=1e-3,
    monitor="pr_auc_score",
    monitor_criterion="max"
)

# 7. Evaluate
print("Evaluating...")
results = trainer.evaluate(test_loader)
print(json.dumps(results, indent=2))
'''
    
    with open(output_path, 'w') as f:
        f.write(template)
    
    print(f"Mortality prediction template saved to {output_path}")


def generate_drug_recommendation_template(output_path: str):
    """Generate a template script for drug recommendation."""
    template = '''#!/usr/bin/env python3
"""
Drug Recommendation using PyHealth
Example: SafeDrug for medication recommendation
"""

from pyhealth.datasets import MIMIC4Dataset
from pyhealth.tasks import drug_recommendation_mimic4_fn
from pyhealth.datasets import split_by_patient, get_dataloader
from pyhealth.models import SafeDrug
from pyhealth.trainer import Trainer

# 1. Load dataset
print("Loading MIMIC-IV dataset...")
dataset = MIMIC4Dataset(
    root="/path/to/mimic4",
    tables=["diagnoses", "procedures", "prescriptions"],
    dev=True
)

# 2. Set task
print("Setting drug recommendation task...")
sample_dataset = dataset.set_task(drug_recommendation_mimic4_fn)

# 3. Split data
train_dataset, val_dataset, test_dataset = split_by_patient(
    sample_dataset, ratios=[0.7, 0.1, 0.2]
)

# 4. Create dataloaders
train_loader = get_dataloader(train_dataset, batch_size=64, shuffle=True)
val_loader = get_dataloader(val_dataset, batch_size=64)
test_loader = get_dataloader(test_dataset, batch_size=64)

# 5. Initialize SafeDrug model
print("Initializing SafeDrug model...")
model = SafeDrug(
    dataset=sample_dataset,
    feature_keys=["conditions", "procedures"],
    label_key="drugs",
    mode="multilabel",
    embedding_dim=64,
    hidden_dim=64
)

# 6. Train
trainer = Trainer(model=model, device="cuda")
trainer.train(
    train_dataloader=train_loader,
    val_dataloader=val_loader,
    epochs=50,
    monitor="jaccard_score"
)

# 7. Evaluate
results = trainer.evaluate(test_loader)
print(results)
'''
    
    with open(output_path, 'w') as f:
        f.write(template)
    
    print(f"Drug recommendation template saved to {output_path}")


def translate_medical_code(code: str, from_system: str, to_system: str) -> str:
    """Translate medical codes between coding systems."""
    try:
        from pyhealth.medical_code_mapper import CrossMap
        
        mapping = CrossMap(from_system, to_system)
        result = mapping.map(code)
        
        return result[0] if result else ""
    except ImportError:
        print("PyHealth not installed")
        return ""
    except Exception as e:
        print(f"Error translating code: {e}")
        return ""


def get_coding_systems() -> List[Dict]:
    """List supported medical coding systems."""
    systems = [
        {"code": "ICD9CM", "name": "ICD-9-CM", "type": "Diagnosis/Procedure"},
        {"code": "ICD10CM", "name": "ICD-10-CM", "type": "Diagnosis"},
        {"code": "ICD10PCS", "name": "ICD-10-PCS", "type": "Procedure"},
        {"code": "CCSCM", "name": "CCS Clinical", "type": "Diagnosis Category"},
        {"code": "CCSPROC", "name": "CCS Procedure", "type": "Procedure Category"},
        {"code": "ATC", "name": "ATC", "type": "Drug Classification"},
        {"code": "NDC", "name": "NDC", "type": "Drug Code"},
        {"code": "RXNORM", "name": "RxNorm", "type": "Drug Normalization"},
    ]
    return systems


def calculate_metrics(y_true: List, y_pred: List, task_type: str = "binary") -> Dict:
    """Calculate evaluation metrics for healthcare predictions."""
    try:
        from sklearn.metrics import (
            accuracy_score, precision_score, recall_score, f1_score,
            roc_auc_score, average_precision_score, 
            jaccard_score, mean_absolute_error
        )
        
        metrics = {}
        
        if task_type == "binary":
            metrics["accuracy"] = accuracy_score(y_true, y_pred > 0.5)
            metrics["precision"] = precision_score(y_true, y_pred > 0.5)
            metrics["recall"] = recall_score(y_true, y_pred > 0.5)
            metrics["f1"] = f1_score(y_true, y_pred > 0.5)
            
            if len(set(y_true)) > 1:
                metrics["auroc"] = roc_auc_score(y_true, y_pred)
                metrics["auprc"] = average_precision_score(y_true, y_pred)
        
        elif task_type == "multilabel":
            metrics["jaccard"] = jaccard_score(y_true, y_pred > 0.5, average="macro")
            metrics["f1_micro"] = f1_score(y_true, y_pred > 0.5, average="micro")
            metrics["f1_macro"] = f1_score(y_true, y_pred > 0.5, average="macro")
        
        return metrics
    except ImportError:
        print("scikit-learn not installed")
        return {}


def create_project_structure(project_name: str):
    """Create a PyHealth project directory structure."""
    import os
    
    dirs = [
        f"{project_name}/data",
        f"{project_name}/models",
        f"{project_name}/scripts",
        f"{project_name}/results",
        f"{project_name}/configs"
    ]
    
    for d in dirs:
        os.makedirs(d, exist_ok=True)
    
    # Create config file
    config = {
        "dataset": {
            "name": "MIMIC4Dataset",
            "root": "./data/mimic4",
            "tables": ["diagnoses", "procedures", "prescriptions"]
        },
        "task": {
            "name": "mortality_prediction",
            "mode": "binary"
        },
        "model": {
            "name": "RETAIN",
            "embedding_dim": 128,
            "hidden_dim": 128
        },
        "training": {
            "batch_size": 64,
            "epochs": 50,
            "learning_rate": 0.001,
            "device": "cuda"
        }
    }
    
    with open(f"{project_name}/configs/config.json", 'w') as f:
        json.dump(config, f, indent=2)
    
    print(f"Project structure created: {project_name}/")


def main():
    parser = argparse.ArgumentParser(description="PyHealth Healthcare AI Tools")
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # List datasets command
    subparsers.add_parser("list-datasets", help="List available datasets")
    
    # List tasks command
    subparsers.add_parser("list-tasks", help="List prediction tasks")
    
    # List models command
    subparsers.add_parser("list-models", help="List available models")
    
    # List coding systems command
    subparsers.add_parser("list-coding-systems", help="List medical coding systems")
    
    # Template commands
    mortality_parser = subparsers.add_parser("template-mortality", help="Generate mortality prediction template")
    mortality_parser.add_argument("-o", "--output", default="mortality_prediction.py")
    
    drug_parser = subparsers.add_parser("template-drug", help="Generate drug recommendation template")
    drug_parser.add_argument("-o", "--output", default="drug_recommendation.py")
    
    # Translate code command
    translate_parser = subparsers.add_parser("translate-code", help="Translate medical code")
    translate_parser.add_argument("code", help="Code to translate")
    translate_parser.add_argument("--from", dest="from_system", required=True, help="Source coding system")
    translate_parser.add_argument("--to", dest="to_system", required=True, help="Target coding system")
    
    # Calculate metrics command
    metrics_parser = subparsers.add_parser("calculate-metrics", help="Calculate evaluation metrics")
    metrics_parser.add_argument("--true", required=True, help="True labels (JSON array)")
    metrics_parser.add_argument("--pred", required=True, help="Predicted probabilities (JSON array)")
    metrics_parser.add_argument("--type", choices=["binary", "multiclass", "multilabel"], default="binary")
    
    # Create project command
    project_parser = subparsers.add_parser("create-project", help="Create project structure")
    project_parser.add_argument("name", help="Project name")
    
    args = parser.parse_args()
    
    if args.command == "list-datasets":
        datasets = list_available_datasets()
        print(json.dumps(datasets, indent=2))
    
    elif args.command == "list-tasks":
        tasks = list_prediction_tasks()
        print(json.dumps(tasks, indent=2))
    
    elif args.command == "list-models":
        models = list_models()
        print(json.dumps(models, indent=2))
    
    elif args.command == "list-coding-systems":
        systems = get_coding_systems()
        print(json.dumps(systems, indent=2))
    
    elif args.command == "template-mortality":
        generate_mortality_prediction_template(args.output)
    
    elif args.command == "template-drug":
        generate_drug_recommendation_template(args.output)
    
    elif args.command == "translate-code":
        result = translate_medical_code(args.code, args.from_system, args.to_system)
        print(f"{args.code} ({args.from_system}) -> {result} ({args.to_system})")
    
    elif args.command == "calculate-metrics":
        y_true = json.loads(args.true)
        y_pred = json.loads(args.pred)
        metrics = calculate_metrics(y_true, y_pred, args.type)
        print(json.dumps(metrics, indent=2))
    
    elif args.command == "create-project":
        create_project_structure(args.name)
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
