import sqlite3

class clique_chunk_summary:
    """A summary of clique gene presence in cliques"""

    def __init__(self, clique_genes, clique_size, startpoint, n_of_cliques):
        """Create a summary object to be stored in the clque database

        genes dictionary with genes as keys and number
        clique_size size of the clique and so the table it is associated with
        startpoint  The starting point of the chunk that is summariazed
        n_of_cliques How many cliques have been evaluated
        """
        self._genes = {}

        self._clique_size = clique_size
        self._startpoint = startpoint
        self._n_of_cliques = n_of_cliques

        self.update_genes(clique_genes)


    def update_genes(self, clique_genes):
        for clique_gene in clique_genes:
            if clique_gene in self._genes:
                self._genes[clique_gene] += 1
            else:
                self._genes[clique_gene] = 1

    def update_n_of_cliques(self):
        self._n_of_cliques += 1

    def get_genes(self):
        return self._genes

    def get_n_of_cliques(self):
        return self._n_of_cliques

    def get_startpoint(self):
        return self._startpoint

    def get_clique_size(self):
        return (self._clique_size)

    def to_db(self, cursor):
        query = ', '.join(["INSERT INTO chunk_summary VALUES (" + \
        str(self.get_startpoint()), str(self.get_startpoint() + \
        self.get_n_of_cliques() -1 ), str(self.get_clique_size()),  \
        str("'" + ' '.join(self.get_genes().keys())) + "');"])
        cursor.execute(query)

gene_presence = {}
unique_genes = set()


def read_temp_file(tempfile, vertex_dictionary, c):
    for line in open(tempfile, "r"):
        insert_into_table(line, c, vertex_dictionary)

def create_summary_table(c):
    query =  '''CREATE TABLE IF NOT EXISTS chunk_summary (
                start INTEGER,
                stop INTEGER,
                clique_size INTEGER,
                genes TEXT);'''

    c.execute(query)

def connect(db_folder, db_file):
    con = sqlite3.connect(db_folder + db_file)
    return con


def insert_into_table(clique, c, vertex_dictionary):
    global unique_genes
    clique_split = clique.rstrip().split(" ")
    genes = replace_index_by_gene(clique_split, vertex_dictionary)
    unique_genes.update(genes)
    table_name = "clique" + str(len(clique_split))
    if table_name in gene_presence:
        if gene_presence[table_name].get_n_of_cliques() == 200000:
            gene_presence[table_name].to_db(c)
            new_clique_summary = clique_chunk_summary(genes, \
            len(genes), \
            (gene_presence[table_name].get_startpoint() +  \
            gene_presence[table_name].get_n_of_cliques()), 1)
            gene_presence[table_name] = new_clique_summary
        else:
            gene_presence[table_name].update_genes(genes)
            gene_presence[table_name].update_n_of_cliques()

    else:
        new_clique_summary = clique_chunk_summary(genes, len(genes), 1, 1)
        gene_presence[table_name] = new_clique_summary

    #If table already exists move on
    if (check_table_exists(table_name, c)):
        pass
    #If not, build the CREATE TABLE query and create the table
    else:
        gene_columns = ""
        for gene_n in range(len(clique_split)):
            gene_columns = gene_columns + " gene" + str((gene_n+1)) + " TEXT, "
        gene_columns = gene_columns[1:-2]

        sql_table = " CREATE TABLE IF NOT EXISTS " + table_name + " (" + gene_columns + ");"
        c.execute(sql_table)
    #Insert the clique into the correct table
    query = "INSERT INTO " + table_name + " VALUES (" + ', '.join(str(x) for x in genes) + ");"
    c.execute(query)

def replace_index_by_gene(clique_split, vertex_dictionary):
    genes = []
    for gene in clique_split:
        gene = str(int(gene))
        genes.append(vertex_dictionary.get(str(int(gene))))

    return(genes)



def check_table_exists(table_name, c):
    query = "SELECT name FROM sqlite_master WHERE type='table' AND name= '" + table_name + "';"
    c.execute(query)

    if (c.fetchone()) !=  None:
        return True
    else:
        return False



def build_clique_db(db_folder, dbname, tempfile, vertex_dictionary_file):
    vertex_dictionary = {}
    for line in open(vertex_dictionary_file, "r"):
        line = line.rstrip().split(' ')
        vertex_dictionary[line[0]] = line[1]
    conn = connect(db_folder, dbname)
    c = conn.cursor()
    create_summary_table(c)
    read_temp_file(tempfile, vertex_dictionary, c)
    for k, v in gene_presence.items():
        v.to_db(c)
    c.execute("CREATE TABLE IF NOT EXISTS unique_genes (gene TEXT);")
    for gene in unique_genes:
        query = "INSERT INTO unique_genes VALUES (" + gene + ");"
        c.execute(query)
    conn.commit()
    c.close()
    conn.close()
