#!/usr/bin/env python
# -*- coding: utf-8 -*-

import mysql.connector
from mysql.connector import Error

from getpass import getpass

nustruDB = mysql.connector.connect(
    host="localhost",
    user=input("Enter username: "),
    password=getpass("Enter password: "),
    database="nustruDB"
)

def create_table():
    """
    Create database connection to nustruDB
    """

    try:   
        sqlCursor = nustruDB.cursor()
        
        # Create a table
        sqlCursor.execute('''CREATE TABLE IF NOT EXISTS nucleotide_protein_seqs
                        (id INT AUTO_INCREMENT PRIMARY KEY,
                        primary_id CHAR(20),
                        gene_name CHAR(80),
                        organism CHAR(200),
                        expression_system CHAR(200),
                        mitochondrial enum('True','False'),
                        protein_sequence MEDIUMTEXT,
                        nucleotide_id CHAR(20), 
                        nucleotide_sequence MEDIUMTEXT);''')
        
        nustruDB.commit()
    
    except Error as e:
        print('Error occurred - ', e)
        
    
create_table()