#!/usr/bin/env python3
"""
Transformers: Hugging Face NLP Examples
========================================
Pre-trained transformer models for NLP, vision, and more.
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np

# Check for transformers availability
try:
    from transformers import (
        pipeline, AutoTokenizer, AutoModel, AutoModelForCausalLM,
        AutoModelForSequenceClassification, AutoModelForQuestionAnswering,
        AutoModelForTokenClassification, Trainer, TrainingArguments
    )
    from transformers import set_seed
    TRANSFORMERS_AVAILABLE = True
except ImportError:
    TRANSFORMERS_AVAILABLE = False
    print("Warning: transformers not installed. Some examples will be skipped.")


# ==============================================================================
# Example 1: Text Classification with Pipelines
# ==============================================================================

def example_text_classification():
    """Sentiment analysis and text classification using pipelines."""
    print("=" * 60)
    print("Example 1: Text Classification (Sentiment Analysis)")
    print("=" * 60)
    
    if not TRANSFORMERS_AVAILABLE:
        print("Transformers not available. Skipping.")
        return None
    
    try:
        # Create sentiment analysis pipeline
        classifier = pipeline(
            "sentiment-analysis",
            model="distilbert-base-uncased-finetuned-sst-2-english"
        )
        
        # Test texts
        texts = [
            "I love this product! It's amazing.",
            "This is the worst experience ever.",
            "The movie was okay, nothing special.",
            "Absolutely fantastic! Highly recommended.",
            "I'm disappointed with the quality."
        ]
        
        print("\nSentiment Analysis Results:")
        print("-" * 50)
        
        for text in texts:
            result = classifier(text)[0]
            print(f"Text: {text[:50]}...")
            print(f"  Label: {result['label']}, Score: {result['score']:.4f}")
            print()
        
        return classifier
        
    except Exception as e:
        print(f"Error: {e}")
        return None


# ==============================================================================
# Example 2: Text Generation
# ==============================================================================

def example_text_generation():
    """Generate text with causal language models."""
    print("\n" + "=" * 60)
    print("Example 2: Text Generation")
    print("=" * 60)
    
    if not TRANSFORMERS_AVAILABLE:
        print("Transformers not available. Skipping.")
        return None
    
    try:
        # Create text generation pipeline
        generator = pipeline(
            "text-generation",
            model="gpt2",
            device=-1  # CPU
        )
        
        # Prompts
        prompts = [
            "The future of artificial intelligence is",
            "In the field of drug discovery, machine learning",
            "Climate change research requires"
        ]
        
        print("\nText Generation Results:")
        print("-" * 50)
        
        for prompt in prompts:
            print(f"\nPrompt: '{prompt}'")
            results = generator(
                prompt,
                max_length=50,
                num_return_sequences=1,
                temperature=0.8,
                do_sample=True
            )
            generated = results[0]['generated_text']
            print(f"Generated: {generated}")
            print()
        
        return generator
        
    except Exception as e:
        print(f"Error: {e}")
        return None


# ==============================================================================
# Example 3: Question Answering
# ==============================================================================

def example_question_answering():
    """Extractive question answering."""
    print("\n" + "=" * 60)
    print("Example 3: Question Answering")
    print("=" * 60)
    
    if not TRANSFORMERS_AVAILABLE:
        print("Transformers not available. Skipping.")
        return None
    
    try:
        # Create QA pipeline
        qa_pipeline = pipeline(
            "question-answering",
            model="distilbert-base-cased-distilled-squad"
        )
        
        # Context and questions
        context = """
        Transformers is a library produced by Hugging Face that provides 
        general-purpose architectures for Natural Language Understanding (NLU) 
        and Natural Language Generation (NLG) with over 32 pretrained models 
        in 100+ languages. It enables fine-tuning of models for various tasks 
        including text classification, question answering, and text generation.
        """
        
        questions = [
            "What is Transformers?",
            "Who produces the Transformers library?",
            "How many pretrained models are available?",
            "What tasks can be performed?"
        ]
        
        print("\nQuestion Answering Results:")
        print("-" * 50)
        
        for question in questions:
            result = qa_pipeline(question=question, context=context)
            print(f"Q: {question}")
            print(f"A: {result['answer']}")
            print(f"  Score: {result['score']:.4f}")
            print()
        
        return qa_pipeline
        
    except Exception as e:
        print(f"Error: {e}")
        return None


# ==============================================================================
# Example 4: Named Entity Recognition (NER)
# ==============================================================================

def example_ner():
    """Named entity recognition."""
    print("\n" + "=" * 60)
    print("Example 4: Named Entity Recognition")
    print("=" * 60)
    
    if not TRANSFORMERS_AVAILABLE:
        print("Transformers not available. Skipping.")
        return None
    
    try:
        # Create NER pipeline
        ner_pipeline = pipeline(
            "ner",
            model="dslim/bert-base-NER",
            aggregation_strategy="simple"
        )
        
        # Test text
        texts = [
            "Apple is looking at buying U.K. startup for $1 billion",
            "Barack Obama was born in Hawaii and was president of the United States.",
            "Transformers was developed by Hugging Face, a company based in New York."
        ]
        
        print("\nNamed Entity Recognition Results:")
        print("-" * 50)
        
        for text in texts:
            print(f"\nText: {text}")
            entities = ner_pipeline(text)
            for ent in entities:
                print(f"  - {ent['word']} ({ent['entity_group']}): {ent['score']:.4f}")
        
        return ner_pipeline
        
    except Exception as e:
        print(f"Error: {e}")
        return None


# ==============================================================================
# Example 5: Feature Extraction (Embeddings)
# ==============================================================================

def example_feature_extraction():
    """Extract embeddings from transformer models."""
    print("\n" + "=" * 60)
    print("Example 5: Feature Extraction (Embeddings)")
    print("=" * 60)
    
    if not TRANSFORMERS_AVAILABLE:
        print("Transformers not available. Skipping.")
        return None
    
    try:
        # Load model and tokenizer
        model_name = "sentence-transformers/all-MiniLM-L6-v2"
        tokenizer = AutoTokenizer.from_pretrained(model_name)
        model = AutoModel.from_pretrained(model_name)
        
        # Texts to encode
        texts = [
            "Machine learning is fascinating.",
            "Deep learning transforms industries.",
            "Neural networks are powerful models."
        ]
        
        print(f"\nModel: {model_name}")
        print("\nComputing embeddings...")
        
        # Tokenize
        inputs = tokenizer(
            texts,
            padding=True,
            truncation=True,
            return_tensors="pt"
        )
        
        # Get embeddings
        import torch
        with torch.no_grad():
            outputs = model(**inputs)
            # Mean pooling
            embeddings = outputs.last_hidden_state.mean(dim=1)
        
        print(f"Embedding shape: {embeddings.shape}")
        
        # Compute similarities
        from sklearn.metrics.pairwise import cosine_similarity
        sim_matrix = cosine_similarity(embeddings.numpy())
        
        print("\nCosine Similarities:")
        print("-" * 50)
        for i in range(len(texts)):
            for j in range(i+1, len(texts)):
                print(f"Text {i+1} <-> Text {j+1}: {sim_matrix[i, j]:.4f}")
        
        return embeddings
        
    except Exception as e:
        print(f"Error: {e}")
        return None


# ==============================================================================
# Example 6: Tokenization
# ==============================================================================

def example_tokenization():
    """Understanding tokenization."""
    print("\n" + "=" * 60)
    print("Example 6: Tokenization")
    print("=" * 60)
    
    if not TRANSFORMERS_AVAILABLE:
        print("Transformers not available. Skipping.")
        return None
    
    try:
        # Load tokenizer
        tokenizer = AutoTokenizer.from_pretrained("bert-base-uncased")
        
        texts = [
            "Hello, how are you?",
            "Transformers are powerful models for NLP.",
            "Tokenization splits text into subword units."
        ]
        
        print(f"\nTokenizer: bert-base-uncased")
        print(f"Vocab size: {tokenizer.vocab_size}")
        
        for text in texts:
            print(f"\nText: '{text}'")
            
            # Tokenize
            tokens = tokenizer.tokenize(text)
            print(f"  Tokens: {tokens}")
            
            # Convert to IDs
            token_ids = tokenizer.convert_tokens_to_ids(tokens)
            print(f"  Token IDs: {token_ids}")
            
            # Encode (with special tokens)
            encoded = tokenizer.encode(text, add_special_tokens=True)
            print(f"  Encoded (with special tokens): {encoded}")
            
            # Decode back
            decoded = tokenizer.decode(encoded)
            print(f"  Decoded: '{decoded}'")
        
        # Padding and truncation
        batch_texts = [
            "Short.",
            "This is a much longer sentence that might need truncation.",
            "Medium length sentence here."
        ]
        
        print("\nBatch encoding with padding:")
        batch_encoded = tokenizer(
            batch_texts,
            padding=True,
            truncation=True,
            max_length=20,
            return_tensors="pt"
        )
        
        print(f"  Input IDs shape: {batch_encoded['input_ids'].shape}")
        print(f"  Attention mask shape: {batch_encoded['attention_mask'].shape}")
        
        return tokenizer
        
    except Exception as e:
        print(f"Error: {e}")
        return None


# ==============================================================================
# Example 7: Summarization
# ==============================================================================

def example_summarization():
    """Text summarization."""
    print("\n" + "=" * 60)
    print("Example 7: Text Summarization")
    print("=" * 60)
    
    if not TRANSFORMERS_AVAILABLE:
        print("Transformers not available. Skipping.")
        return None
    
    try:
        # Create summarization pipeline
        summarizer = pipeline(
            "summarization",
            model="facebook/bart-large-cnn"
        )
        
        # Long text to summarize
        text = """
        Transformers have revolutionized natural language processing since the 
        introduction of the "Attention Is All You Need" paper in 2017. The 
        architecture, based solely on attention mechanisms, has enabled unprecedented 
        scale in language models, leading to breakthroughs in understanding and 
        generating human language. Models like BERT, GPT, and T5 have set new 
        benchmarks across a wide range of NLP tasks, from sentiment analysis to 
        machine translation. The key innovation of the transformer is the 
        self-attention mechanism, which allows the model to weigh the importance 
        of different words in a sentence, regardless of their position. This has 
        overcome the limitations of previous sequential models like RNNs and LSTMs, 
        which struggled with long-range dependencies. Today, transformers are not 
        only used for text but have been adapted for images, audio, and multimodal 
        tasks, demonstrating their versatility as a general-purpose architecture 
        for deep learning.
        """
        
        print(f"\nOriginal text length: {len(text)} characters")
        
        # Generate summary
        summary = summarizer(
            text,
            max_length=100,
            min_length=30,
            do_sample=False
        )
        
        print(f"\nSummary ({len(summary[0]['summary_text'])} chars):")
        print("-" * 50)
        print(summary[0]['summary_text'])
        
        return summarizer
        
    except Exception as e:
        print(f"Error: {e}")
        return None


# ==============================================================================
# Example 8: Translation
# ==============================================================================

def example_translation():
    """Machine translation."""
    print("\n" + "=" * 60)
    print("Example 8: Machine Translation")
    print("=" * 60)
    
    if not TRANSFORMERS_AVAILABLE:
        print("Transformers not available. Skipping.")
        return None
    
    try:
        # Create translation pipeline
        translator = pipeline(
            "translation",
            model="Helsinki-NLP/opus-mt-en-de"
        )
        
        # English texts
        texts = [
            "Hello, how are you today?",
            "Machine learning is transforming the world.",
            "I love learning new languages."
        ]
        
        print("\nEnglish -> German Translation:")
        print("-" * 50)
        
        for text in texts:
            result = translator(text, max_length=100)
            print(f"EN: {text}")
            print(f"DE: {result[0]['translation_text']}")
            print()
        
        return translator
        
    except Exception as e:
        print(f"Error: {e}")
        return None


# ==============================================================================
# Example 9: Fill Mask (Masked Language Modeling)
# ==============================================================================

def example_fill_mask():
    """Masked language modeling."""
    print("\n" + "=" * 60)
    print("Example 9: Fill Mask (Masked Language Modeling)")
    print("=" * 60)
    
    if not TRANSFORMERS_AVAILABLE:
        print("Transformers not available. Skipping.")
        return None
    
    try:
        # Create fill-mask pipeline
        unmasker = pipeline(
            "fill-mask",
            model="bert-base-uncased"
        )
        
        # Masked texts
        texts = [
            "The capital of France is [MASK].",
            "Machine learning is a subset of [MASK].",
            "The [MASK] is the largest planet in our solar system."
        ]
        
        print("\nFill Mask Predictions:")
        print("-" * 50)
        
        for text in texts:
            print(f"\nText: '{text}'")
            predictions = unmasker(text, top_k=3)
            for pred in predictions:
                print(f"  {pred['token_str']}: {pred['score']:.4f}")
        
        return unmasker
        
    except Exception as e:
        print(f"Error: {e}")
        return None


# ==============================================================================
# Main Execution
# ==============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("Transformers: Hugging Face NLP Examples")
    print("=" * 70)
    
    if not TRANSFORMERS_AVAILABLE:
        print("\nNote: Transformers library is not installed.")
        print("To install: pip install transformers torch")
        print("Examples will run in demonstration mode.\n")
    
    try:
        example_text_classification()
        example_text_generation()
        example_question_answering()
        example_ner()
        example_feature_extraction()
        example_tokenization()
        example_summarization()
        example_translation()
        example_fill_mask()
        
        print("\n" + "=" * 70)
        print("All examples completed!")
        print("=" * 70)
        
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
