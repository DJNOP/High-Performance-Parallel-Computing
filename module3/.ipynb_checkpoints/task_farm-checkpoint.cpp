/*
  Assignment: Make an MPI task farm. A "task" is a randomly generated integer.
  To "execute" a task, the worker sleeps for the given number of milliseconds.
  The result of a task should be send back from the worker to the master. It
  contains the rank of the worker
*/

#include <iostream>
#include <random>
#include <chrono>
#include <thread>
#include <array>

// To run an MPI program we always need to include the MPI headers
#include <mpi.h>

const int NTASKS=5000;  // number of tasks
const int RANDOM_SEED=1234;

void master (int nworker) {
    std::array<int, NTASKS> task, result;

    // set up a random number generator
    std::random_device rd;
    //std::default_random_engine engine(rd());
    std::default_random_engine engine;
    engine.seed(RANDOM_SEED);
    // make a distribution of random integers in the interval [0:30]
    std::uniform_int_distribution<int> distribution(0, 30);

    for (int& t : task) {
        t = distribution(engine);   // set up some "tasks"
    }

    result.fill(-1);

    const int TAG_WORK   = 1;  // master -> worker: {task_id, task_ms}
    const int TAG_STOP   = 2;  // master -> worker: {-1, 0}
    const int TAG_RESULT = 3;  // worker -> master: {task_id}

    int next_task = 0;  // next task index to hand out
    int done      = 0;  // number of completed tasks received

    // initial dispatch: at most one task per worker
    for (int worker = 1; worker <= nworker; ++worker) {
        if (next_task < NTASKS) {
            int msg[2] = { next_task, task[next_task] };
            MPI_Send(msg, 2, MPI_INT, worker, TAG_WORK, MPI_COMM_WORLD);
            ++next_task;
        } else {
            int msg[2] = { -1, 0 };
            MPI_Send(msg, 2, MPI_INT, worker, TAG_STOP, MPI_COMM_WORLD);
        }
    }

    // main loop: take result from whoever finishes, then give that worker new work (or stop)
    while (done < NTASKS) {
        MPI_Status status;
        int task_id = -1;

        MPI_Recv(&task_id, 1, MPI_INT, MPI_ANY_SOURCE, TAG_RESULT,
                 MPI_COMM_WORLD, &status);

        const int worker = status.MPI_SOURCE;

        if (0 <= task_id && task_id < NTASKS) {
            result[task_id] = worker;  // store who did this task
            ++done;
        }

        if (next_task < NTASKS) {
            int msg[2] = { next_task, task[next_task] };
            MPI_Send(msg, 2, MPI_INT, worker, TAG_WORK, MPI_COMM_WORLD);
            ++next_task;
        } else {
            int msg[2] = { -1, 0 };
            MPI_Send(msg, 2, MPI_INT, worker, TAG_STOP, MPI_COMM_WORLD);
        }
    }

    // Print out a status on how many tasks were completed by each worker
    for (int worker=1; worker<=nworker; worker++) {
        int tasksdone = 0; int workdone = 0;
        for (int itask=0; itask<NTASKS; itask++)
        if (result[itask]==worker) {
            tasksdone++;
            workdone += task[itask];
        }
        std::cout << "Master: Worker " << worker << " solved " << tasksdone <<
                    " tasks\n";
    }
}

// call this function to complete the task. It sleeps for task milliseconds
void task_function(int task) {
    std::this_thread::sleep_for(std::chrono::milliseconds(task));
}

void worker (int rank) {
    const int TAG_WORK   = 1;  // master -> worker: {task_id, task_ms}
    const int TAG_STOP   = 2;  // master -> worker: {-1, 0}
    const int TAG_RESULT = 3;  // worker -> master: {task_id}

    while (true) {
        int msg[2] = { -1, 0 };
        MPI_Status status;

        // Wait for either a new task or a stop signal from the master
        MPI_Recv(msg, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        const int task_id = msg[0];
        const int task_ms = msg[1];

        // Stop condition: master sends {-1, 0} with TAG_STOP
        // (also safe if you just check task_id == -1)
        if (status.MPI_TAG == TAG_STOP || task_id == -1) {
            break;
        }

        // Do the work (sleep for task_ms milliseconds)
        task_function(task_ms);

        // Send back which task was completed
        MPI_Send(&task_id, 1, MPI_INT, 0, TAG_RESULT, MPI_COMM_WORLD);
    }
}

int main(int argc, char *argv[]) {
    int nrank, rank;

    MPI_Init(&argc, &argv);                // set up MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nrank); // get the total number of ranks
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // get the rank of this process

    if (rank == 0)       // rank 0 is the master
        master(nrank-1); // there is nrank-1 worker processes
    else                 // ranks in [1:nrank] are workers
        worker(rank);

    MPI_Finalize();      // shutdown MPI
}
