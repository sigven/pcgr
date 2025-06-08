import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from api.validators import VariantValidator
from api.database import get_session

def test_validation():
    session = next(get_session())
    try:
        validator = VariantValidator(session)
        
        test_cases = ["MSI-H", "TMB:15.2", "HRD:45.0", "BRAF_expression:2.5"]
        for case in test_cases:
            result = validator.validate_batch([case], 'MEL')
            print(f'Testing {case}:')
            print(f'  Errors: {result.summary["errors"]}')
            if result.errors:
                print(f'  Error: {result.errors[0].message}')
            print(f'  Is tumor marker: {validator._is_tumor_marker(case)}')
            print()
        
        result = validator.validate_batch(['NM_004333.4:c.1799T>A', 'chr7:140753336A>T', 'MSI-H'], 'MEL')
        print('Original validation test result:')
        print(f'Total: {result.summary["total_variants"]}')
        print(f'Parsed: {result.summary["parsed"]}')
        print(f'Errors: {result.summary["errors"]}')
        if result.errors:
            for error in result.errors:
                print(f'Error line {error.line}: {error.message}')
    finally:
        session.close()

if __name__ == "__main__":
    test_validation()
