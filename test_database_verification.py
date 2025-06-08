import sqlite3
import os

def verify_database():
    """Verify database population after ETL pipeline."""
    db_path = 'dgg_rules_demo.db'
    if not os.path.exists(db_path):
        print(f"Database {db_path} not found")
        return
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    tables = [
        'genegeneralcomment',
        'variantdxinterpretation', 
        'pertinentnegative',
        'variantgeneralcomment',
        'drugassociation'
    ]
    
    for table in tables:
        try:
            cursor.execute(f'SELECT COUNT(*) FROM {table}')
            count = cursor.fetchone()[0]
            print(f'{table}: {count} entries')
        except Exception as e:
            print(f'{table}: Error - {e}')
    
    conn.close()

if __name__ == "__main__":
    verify_database()
