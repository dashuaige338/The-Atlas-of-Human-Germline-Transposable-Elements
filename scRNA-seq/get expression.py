import psycopg2
import pandas as pd
import numpy as np
import json


def get_expression(gene_name):
    conn = psycopg2.connect("dbname = hpgcdb user = postgres password = Hrbl79eY host = 10.106.125.216 port = 5432")
    c = conn.cursor()
    c.execute("select expression from expression where name = '{}'".format(gene_name))
    results=c.fetchall()
    print(results)
    result = results[0][0].tobytes()
    result = json.loads(result)
    result = pd.DataFrame(result,columns=['expression'])
    return result

result = get_expression('ZNF146')
print(result)
