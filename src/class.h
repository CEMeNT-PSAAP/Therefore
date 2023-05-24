class cell{
    public:
        int cell_id;
        float x_left;
        float xsec_scatter;
        float xsec_total;
        float dx;
        float v;
        float dt;
};

class problem_space{
    public:
        // physical space
        int N_cells; 
        int N_angles; 
        int N_time;
        float Length;
        float dt;
        float dx;
        float t_max;
        float material_source;
        float velocity;

        // comptuational
        int hardware_precision;
        float convergence_tolarance;
};

class boundary_condition{
    public:
        char side;
        int cell_id;
        int type;
        float magnitude;
};
