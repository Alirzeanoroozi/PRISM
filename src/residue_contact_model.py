"""Optional residue-contact model for the Stage 2 PRISM experiment.

PyTorch is imported lazily so the production pipeline and CPU smoke tests do
not acquire a deep-learning dependency. This model consumes one candidate's
residue features and a contact edge list and predicts native-like probability
and iRMSD. It is intentionally not wired into the default pipeline yet.
"""

from __future__ import annotations


def build_model(node_features: int, hidden_features: int = 64):
    """Build a small contact-message-passing model when PyTorch is available."""
    try:
        import torch
        from torch import nn
    except ImportError as exc:  # pragma: no cover - depends on optional env
        raise RuntimeError(
            "Stage 2 requires an environment with PyTorch; production PRISM "
            "does not install it automatically"
        ) from exc

    if node_features <= 0 or hidden_features <= 0:
        raise ValueError("node_features and hidden_features must be positive")

    class ContactRanker(nn.Module):
        def __init__(self):
            super().__init__()
            self.encoder = nn.Sequential(
                nn.Linear(node_features, hidden_features),
                nn.ReLU(),
                nn.Linear(hidden_features, hidden_features),
                nn.ReLU(),
            )
            self.classifier = nn.Linear(hidden_features, 1)
            self.irmsd = nn.Linear(hidden_features, 1)

        def forward(self, node_features_tensor, edge_index):
            if edge_index.ndim != 2 or edge_index.shape[0] != 2:
                raise ValueError("edge_index must have shape [2, edge_count]")
            encoded = self.encoder(node_features_tensor)
            source, target = edge_index
            messages = encoded[source]
            aggregated = encoded.clone()
            aggregated.index_add_(0, target, messages)
            degree = torch.bincount(target, minlength=encoded.shape[0]).clamp_min(1)
            aggregated = aggregated / degree.to(encoded.dtype).unsqueeze(1)
            pooled = aggregated.mean(dim=0)
            return {
                "native_like_logit": self.classifier(pooled).squeeze(-1),
                "irmsd": self.irmsd(pooled).squeeze(-1),
            }

    return ContactRanker()
