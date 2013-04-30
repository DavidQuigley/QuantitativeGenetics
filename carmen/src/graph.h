enum {
	NODE_GENE = 0,
	NODE_LOCUS = 1
};

typedef boost::unordered_map<int, boost::unordered_map<int, double>* >::iterator graph_itr;
typedef std::vector<std::vector<int>*>* clique_list;

class Graph{
public:
    Graph();
    ~Graph();

	void add_edge( int node_a, int node_b);
    void add_edge( int node_a, int node_b, double weight );
    void add_edge( int node_a, int node_b, double weight, int node_a_type, int node_b_type );
    void add_qtls( std::vector<QTL*>& qtls );
    void clear();
    void clear_redundant( std::vector<std::vector<int>*>* results, std::vector<int>* target );
    void cliques( int min_k, int max_k, std::vector<std::vector<int>*>* results, bool verbose);
    void copy(Graph* g);
    int degree( int node_a );
    void diff( std::vector<int>* A, std::vector<int>* B, std::vector<int>* C); // C = in B but not A
    double edge_weight(int node_a, int node_b);
    void edges( std::vector<std::vector<int>*>* edges);
    void evaluate_quad( boost::unordered_map<int, int>& status, int c, std::vector<int>& neighbors );
    bool has_node( int node_a );
    bool has_edge( int node_a, int node_b, double& weight ); // if present, weight written to parameter
    bool has_locus_neighbor( int node_a );
    bool has_node_label(int node_a);
    void hubs(std::vector<int>& h, int percentile);
    int get_node_index_by_label(std::string label);
    std::string get_node_label_by_index(int node_a);
    bool is_neighbor( int node_a, int node_b );
    void neighbors( int node_a, std::vector<int>& ids);
    void neighbors( int node_a, std::vector<int>& ids, std::vector<double>& weights );
    int node_type(int node_a);
    void nodes( std::vector<int>& n );
    int order(); // number of nodes;
    int n_diff( std::vector<int>* A, std::vector<int>* B);
    void print();
    void print_nodes();
    void remove_edge( int node_a, int node_b );
    void read( std::string filename );
	void read( std::string filename, double min_abs_weight );
	void set_node_label(int node_a, std::string label);
	void set_node_type(int node_a, int node_type);
	void set_verbose( bool verbose );
    void semi_cliques(double completeness, std::vector<std::vector<int>*>* results, bool verbose);
    void shared_neighbors(std::vector<int>* C, std::vector<int>& N);
    void subgraph( std::vector<int>& ids );
    void subgraph( Graph& g, std::vector<int>& ids );  
    void reduce_to_triangles( Graph& g );
    void reduce_to_quads( Graph& g );
	void remove_vectors_by_size( std::vector<std::vector<int>*>* V, int min_size, bool verbose ); // called by BK_maximal_cliques implementation
    void write( std::string fn_out, boost::unordered_map<int, std::string>& translate );
    Graph* clique_graph; // stores output of BK_maximal_cliques

private:
    bool verbose;
	boost::unordered_map<int, boost::unordered_map<int, double>* > adj; // adjacency stored here
	boost::unordered_map<int, int> node_types;
	boost::unordered_map<int, std::string> idx2str;
	boost::unordered_map<std::string, int> str2idx;
};


class Node{
public:
	Node(int id);
	Node* add_child(int id);
	Node* has_child(int id);
	bool has_path(std::vector<int>* path);
	void add_path(std::vector<int>* path);
	void get_kids(std::vector<int>* c);
	void print_kids();
	static const int ROOT = 0;
	~Node();
private:
	int id;
	boost::unordered_map<int, Node*> children;
};

