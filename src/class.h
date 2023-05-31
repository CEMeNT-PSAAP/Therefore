class cell{
    public:
        int cell_id;
        double x_left;
        std::vector<double> xsec_scatter;
        std::vector<double> xsec_total;
        std::vector<double> v;
        double dx;
        double dt;
};

class problem_space{
    public:
        // physical space
        int N_cells; 
        int N_angles; 
        int N_time;
        double Length;
        double dt;
        double dx;
        double t_max;
        double material_source;
        double velocity;

        // comptuational
        int hardware_precision;
        double convergence_tolarance;
};

class boundary_condition{
    public:
        char side;
        int cell_id;
        int type;
        double magnitude;
};
