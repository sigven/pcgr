import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.kb_models import OncoTreeMapping
from api.database import get_session
from datetime import datetime


def populate_oncotree_codes():
    """Populate OncoTree codes from demo UI list."""
    tumor_types = [
        "MEL", "LUAD", "BRCA", "CRC", "PRAD", "BLCA", "HNSC", "KIRC", 
        "LIHC", "THCA", "UCEC", "KIRP", "SARC", "LAML", "PAAD", "GBM",
        "OV", "SKCM", "CESC", "TGCT", "UCS", "CHOL", "THYM", "ACC",
        "MESO", "UVM", "DLBC", "KICH", "READ", "LGG", "PCPG", "STAD"
    ]
    
    session = next(get_session())
    try:
        for code in tumor_types:
            existing = session.query(OncoTreeMapping).filter(
                OncoTreeMapping.oncotree_code == code
            ).first()
            
            if not existing:
                mapping = OncoTreeMapping(
                    oncotree_code=code,
                    disease_name=f"Disease for {code}",
                    created_at=datetime.utcnow(),
                    updated_at=datetime.utcnow(),
                    created_by="system"
                )
                session.add(mapping)
        
        session.commit()
        print(f"Populated {len(tumor_types)} OncoTree codes")
    finally:
        session.close()


if __name__ == "__main__":
    populate_oncotree_codes()
