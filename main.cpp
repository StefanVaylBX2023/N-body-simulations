#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <thread>
#include <mutex>
#include <chrono>
#include <SFML/Graphics.hpp>
#include <queue>
#include <functional>
#include <condition_variable>

const double G = 6.67430e-11;  // Gravitational constant
std::mutex mutex;  



//#######################################################################################
//                       Body Class Implementation
//#######################################################################################

struct Body {
    double mass;
    double x, y;  // Position
    double vx, vy;  // Velocity
    double fx, fy;  // Force

    Body(double m, double xpos, double ypos, double xvel, double yvel)
        : mass(m), x(xpos), y(ypos), vx(xvel), vy(yvel), fx(0.0), fy(0.0) {}

    void updatePosition(double dt) {
        x += vx * dt;
        y += vy * dt;
    }
};



//#######################################################################################
//                       Direct Simulation
//#######################################################################################



// Function to update velocity of a body based on the computed forces
void updateVelocity(Body& body, double dt) {
    double ax = body.fx / body.mass;
    double ay = body.fy / body.mass;

    body.vx += ax * dt;
    body.vy += ay * dt;
}

// Function to compute the gravitational force between two bodies
void computeForce(Body& body1, Body& body2) {
    double dx = body2.x - body1.x;
    double dy = body2.y - body1.y;
    double distSq = dx * dx + dy * dy;
    double dist = std::sqrt(distSq);
    double force = G * body1.mass * body2.mass / (distSq * dist);

    double fx = force * dx / dist;
    double fy = force * dy / dist;

    std::lock_guard<std::mutex> lock(mutex);  // Lock the mutex to protect shared data access
    body1.fx += fx;
    body1.fy += fy;
    body2.fx -= fx;
    body2.fy -= fy;
}

// Function to simulate the dynamics of the bodies for a given range of indices
void simulateRange(std::vector<Body>& bodies, double dt, int start, int end, std::vector<std::vector<Body>>& iterations) {
    for (int i = start; i < end; ++i) {
        for (int j = i + 1; j < bodies.size(); ++j) { //avoid double computation
            computeForce(bodies[i], bodies[j]);
        }
        updateVelocity(bodies[i], dt);
    }

    

    iterations.push_back(bodies);  // Store the iteration
}



// Function to simulate the dynamics of the bodies using multiple threads
void simulate(std::vector<Body>& bodies, double dt, int numThreads, std::vector<std::vector<Body>>& iterations) {
    int n = bodies.size();
    int chunkSize = (n + numThreads - 1) / numThreads;

    std::vector<std::thread> threads;

    for (int t = 0; t < numThreads; ++t) {
        int start = t * chunkSize;
        int end = std::min((t + 1) * chunkSize, n);
        threads.emplace_back(simulateRange, std::ref(bodies), dt, start, end, std::ref(iterations));
    }

    for (auto& thread : threads) {
        thread.join();
    }

    for (int i = 0; i < n; ++i) {
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
        //uncomment to output forces
        //std::cout<< " Force of body " << i << " = " << bodies[i].fx << ", " << bodies[i].fy << ")" << std::endl;
        bodies[i].fx = 0.0;
        bodies[i].fy = 0.0;
    }
}


//#######################################################################################
//                       SFML Visualization
//#######################################################################################

// Function to draw the bodies on the SFML window
void drawBodies(const std::vector<Body>& bodies, sf::RenderWindow& window) {
    window.clear();

    for (const auto& body : bodies) {
        sf::CircleShape circle(body.mass/8000);
        //sf::CircleShape circle(10);
        circle.setPosition(body.x * window.getSize().x, body.y * window.getSize().y);
        window.draw(circle);
    }

    window.display();
}

//#######################################################################################

//                       Auxiliary Functions

//#######################################################################################

void outputVelocities(const std::vector<Body>& bodies) {
    //
    for (int i = 0; i < bodies.size(); ++i) {
        std::cout << "Velocity (vx, vy) of body "<< i << " = " << bodies[i].vx << ", " << bodies[i].vy << ")" << std::endl;
    }
}

void outputPositions(const std::vector<Body>& bodies) {
    for (int i = 0; i < bodies.size(); ++i) {
        std::cout << "Position (x, y) of body "<< i << " = " << bodies[i].x << ", " << bodies[i].y << ")" << std::endl;
    }
}

void outputForces(const std::vector<Body>& bodies) {
    for (int i = 0; i < bodies.size(); ++i) {
        std::cout << "Force (fx, fy) of body "<< i << " = " << bodies[i].fx << ", " << bodies[i].fy << ")" << std::endl;
    }
}



//#######################################################################################

//                       QuadTree Class Implementation

//#######################################################################################
struct QuadTree {
    double x, y;           // Coordinates of the center of the region
    double width, height;  // Dimensions of the region
    QuadTree* nw;          // Pointer to the northwest quadrant
    QuadTree* ne;          // Pointer to the northeast quadrant
    QuadTree* sw;          // Pointer to the southwest quadrant
    QuadTree* se;          // Pointer to the southeast quadrant
    Body* body;            // Pointer to the body contained in the region
    bool hasChildren;      // Flag indicating if the region has children

    QuadTree(double xpos, double ypos, double w, double h)
        : x(xpos), y(ypos), width(w), height(h), nw(nullptr), ne(nullptr), sw(nullptr), se(nullptr), body(nullptr),
          hasChildren(false) {}

    ~QuadTree() {
        delete nw;
        delete ne;
        delete sw;
        delete se;
    }

    void insert(Body* b) {
        //std::cout << "Inserting body at (" << b->x << ", " << b->y << ")" << std::endl;
        if (!hasChildren && body == nullptr) {
            // If the region is empty, store the body in it
            body = b;
        } else {
            if (!hasChildren) {
                // If the region doesn't have children, create them
                subdivide();
            }

            // Insert the body into the appropriate child region
            if (b->x < x && b->y < y) {
                nw->insert(b);
            } else if (b->x >= x && b->y < y) {
                ne->insert(b);
            } else if (b->x < x && b->y >= y) {
                sw->insert(b);
            } else {
                se->insert(b);
            }
        }
    }

    void subdivide() {
        double childWidth = width / 2;
        double childHeight = height / 2;
        //std::cout << "Subdividing region at (" << x << ", " << y << ")" << std::endl;
        nw = new QuadTree(x - childWidth / 2, y - childHeight / 2, childWidth, childHeight);
        ne = new QuadTree(x + childWidth / 2, y - childHeight / 2, childWidth, childHeight);
        sw = new QuadTree(x - childWidth / 2, y + childHeight / 2, childWidth, childHeight);
        se = new QuadTree(x + childWidth / 2, y + childHeight / 2, childWidth, childHeight);

        hasChildren = true;
    }

    void computeDirectForce(Body* body) {
        double dx = this->x - body->x;
        double dy = this->y - body->y;
        double distSq = dx * dx + dy * dy + 0.001;  // Add small value to avoid division by zero
        double force = G * body->mass * this->width * this->height / distSq;
        double angle = std::atan2(dy, dx);
        body->vx += force * std::cos(angle) / body->mass;
        body->vy += force * std::sin(angle) / body->mass;
        //std::cout << "Force on body at (" << body->x << ", " << body->y << ") = (" << body->fx << ", " << body->fy << ")" << std::endl;
        if (this->body != nullptr && this->body != body) {
            double otherForce = G * this->body->mass * this->width * this->height / distSq;
            this->body->vx -= otherForce * std::cos(angle) / this->body->mass;
            this->body->vy -= otherForce * std::sin(angle) / this->body->mass;
        }
    }

    void computeApproximateForce(Body* body) {
        double dx = this->x - body->x;
        double dy = this->y - body->y;
        double dist = std::sqrt(dx * dx + dy * dy + 0.001);  // Add small value to avoid division by zero

        if (this->width / dist < 0.5) {
            computeDirectForce(body);
        } else {
            if (hasChildren) {
                nw->computeApproximateForce(body);
                ne->computeApproximateForce(body);
                sw->computeApproximateForce(body);
                se->computeApproximateForce(body);
            } else if (this->body != nullptr && this->body != body) {
                computeDirectForce(body);
            }
        }
    }

    void computeForces(Body* body) {

        if (hasChildren) {    
            nw->computeApproximateForce(body);
            ne->computeApproximateForce(body);
            sw->computeApproximateForce(body);
            se->computeApproximateForce(body);
        } else if (this->body != nullptr && this->body != body) {
            computeDirectForce(body);
        }
    }
};


class ThreadPool {
public:
    ThreadPool(int numThreads) : stop(false) {
        for (int i = 0; i < numThreads; ++i) {
            workers.emplace_back([this]() {
                while (true) {
                    std::function<void()> task;

                    {
                        std::unique_lock<std::mutex> lock(queueMutex);
                        condition.wait(lock, [this]() { return stop || !tasks.empty(); });

                        if (stop && tasks.empty()) {
                            return;
                        }

                        task = std::move(tasks.front());
                        tasks.pop();
                    }

                    task();
                }
            });
        }
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            stop = true;
        }

        condition.notify_all();

        for (std::thread& worker : workers) {
            worker.join();
        }
    }

    template<class F, class... Args>
    void enqueue(F&& f, Args&&... args) {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            tasks.emplace([f, args...]() { f(args...); });
        }

        condition.notify_one();
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queueMutex;
    std::condition_variable condition;
    bool stop;
};

void computeForces(QuadTree* tree, Body* body) {
    tree->computeForces(body);
}

void updatePosition(Body* body, double dt) {
    body->updatePosition(dt);
}

void simulateTree(std::vector<Body>& bodies, double dt, int numThreads, QuadTree* tree, bool par) {
    if(par == true){
    ThreadPool threadPool(numThreads);

    // Parallelize force computation
    for (auto& body : bodies) {
        threadPool.enqueue(computeForces, tree, &body);
    }

    // Wait for all force computation tasks to finish

    // Parallelize position updates
    for (auto& body : bodies) {
        threadPool.enqueue(updatePosition, &body, dt);
    }
    }else {
    for (auto& body : bodies) {
        tree->computeForces(&body);
    }
   
    // Update positions and velocities
    for (auto& body : bodies) {
        body.updatePosition(dt);
    }
    }
}



int main(int argc, char** argv) {
    std::srand(std::time(0));  // Seed random number generator
    int numSteps = 1000;
    int n = 1000;  // Number of bodies
    double dt = 0.01;  // Time step
    int numThreads = 2;  // Number of threads to use
    bool useBarnesHut = false;  // Flag indicating whether to use Barnes-Hut algorithm
    bool par = false;
    bool naive = false;
    if (argc > 1) {
        std::string algorithm = argv[1];
        if (algorithm == "barnes-hutt-par") {
            useBarnesHut = true;
            par = true;
        }else if (algorithm == "barnes-hutt") {
            useBarnesHut = true;
        }else if (algorithm == "naive") {
            naive = true;
        }
        
    }
    if (argc > 2) {
        n = std::stoi(argv[2]);
    }
    if (argc > 3) {
        numSteps = std::stoi(argv[3]);
    }
    
    if (argc <= 1){
        naive = true;
    }
    // Initialize bodies with random positions, velocities, and masses
    std::vector<Body> bodies;
    bodies.reserve(n);
    for (int i = 0; i < n; ++i) {
        double mass = std::rand() / double(RAND_MAX)*100000;  // Random mass between 0 and 1
        double xpos = (std::rand() / double(RAND_MAX)) * 2 - 1;  // Random x position between -1 and 1
        double ypos = (std::rand() / double(RAND_MAX)) * 2 - 1;  // Random y position between 0 and 1
        double xvel = ((std::rand() /  double(RAND_MAX)) * 2  - 1)/10;  // Random x velocity between 0 and 1
        double yvel = ((std::rand()/  double(RAND_MAX)) * 2  - 1)/10;  // Random y velocity between 0 and 1
        bodies.emplace_back(mass, xpos, ypos, xvel, yvel);
    }

    // Create SFML window
    sf::RenderWindow window(sf::VideoMode(2000, 2000), "N-Body Simulation");
    window.setFramerateLimit(60);

    // Perform simulation for a specified number of steps
    std::vector<std::vector<Body>> iterations;  // Store iterations

    auto start = std::chrono::high_resolution_clock::now();  // Start timer

    if (useBarnesHut) {
        // Perform Barnes-Hut simulation
        QuadTree tree(-1.0, -1.0, 2.0, 2.0);  // Create quadtree covering the entire simulation area
        for (int i = 0; i < n; ++i) {
            tree.insert(&bodies[i]);
        }
        //std::cout << "Inserting bodies into quadtree..." << std::endl;
        for (int step = 0; step < numSteps; ++step) {
            simulateTree(bodies, dt, numThreads, &tree, par);
            iterations.push_back(bodies);
            //std::cout<< bodies[22].vx << ", " << bodies[22].vy << ")" << std::endl;
        }
    } else if (naive) {
        // Perform direct simulation
        std::cout<< "Naive" << std::endl;
        for (int step = 0; step < numSteps; ++step) {
            simulate(bodies, dt, numThreads, iterations);
            //std::cout<< bodies[22].vx << ", " << bodies[22].vy << ")" << std::endl;
        }

    }

    auto end = std::chrono::high_resolution_clock::now();  // End timer
    std::chrono::duration<double> elapsed = end - start;  // Compute elapsed time
    std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;

    // Visualize stored iterations
    for (const auto& iteration : iterations) {
        // Update SFML window
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        drawBodies(iteration, window);
    }

    return 0;

}