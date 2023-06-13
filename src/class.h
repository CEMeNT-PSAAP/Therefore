using namespace std;

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

        // computational
        int hardware_precision;
        double convergence_tolerance;
        bool initialize_from_previous;

        int boundary_condition_left;
        int boundary_condition_right;
        
        double boundary_condition(){
            return(0.0);
        }
};

class boundary_condition{
    public:
        char side;
        int cell_id;
        int type;
        double magnitude;
};


class ts_solutions{
    public:
        std::vector<double> aflux;
        double time;
        double spectral_radius;
        double final_error;
        int number_iteration;
        int N_step;
};

class WholeProblem{
    public:
        vector<cell> cells;
        problem_space ps;
        vector<ts_solutions> solutions;

        WholeProblem(vector<cell> c, problem_space p, vector<ts_solutions> sol){
            cells = c;
            ps = p;
            solutions = sol;
        }

        void PublishWholeProblem(){
            
        }


};