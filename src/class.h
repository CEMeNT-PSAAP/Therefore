class cell{
    public:
        int cell_id;
        double x_left;
        std::vector<double> xsec_scatter;
        std::vector<double> xsec_total;
        std::vector<double> v;
        std::vector<double> Q;
        double dx;
        double dt;
};

class problem_space{
    public:
        // physical space
        int N_cells; 
        int N_angles; 
        int N_time;
        int N_groups;
        double Length;
        double dt;
        double dx;
        double t_max;
        double material_source;
        double velocity;

        double *weights;
        double *angles;

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


class ts_soultion{
    public:
        std::vector<double> aflux;
        std::vector<double> aflux_h;
        std::vector<double> sflux;
        double time;
        double specral_radius;
        int number_itteration;
        int N_step;

};