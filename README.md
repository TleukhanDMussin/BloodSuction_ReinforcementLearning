
# 3D Blood Suction Simulator with Reinforcement Learning

This project demonstrates the integration of a 3D blood suction simulator with a reinforcement learning (RL) agent using Pybind11, Gym, and Stable Baselines3.

## Directory Structure
```
BloodSuctionRL/
│
├── simulator/         # C++ simulation and Pybind11 bindings
├── env/               # Python Gym environment
├── rl/                # RL agent and training
├── assets/            # Original files and assets
└── README.md          # Project README with setup instructions
```

## Requirements
- C++ Compiler (with CMake support)
- Python 3.x
- Libraries: pybind11, gym, numpy, stable-baselines3, matplotlib

## Setup Instructions

### 1. Build the C++ Module
- Go to the `simulator` directory:
  ```bash
  cd simulator
  mkdir build
  cd build
  cmake ..
  make
  ```
- This will create a Python module named `simulator` in the `build` directory.

### 2. Install Python Dependencies
```bash
pip install -r rl/requirements.txt
```

### 3. Run the RL Training
- From the project root, run:
  ```bash
  python rl/train_rl_agent.py
  ```

### 4. Test the Trained Model
- The script will automatically run a test after training and display the simulation.

## Description
This project models the dynamics of blood particles during a surgical suction process and trains an RL agent to control the suction tool to clear all particles from the container.

### Contributions
Feel free to modify, improve, or contribute to the project!

### License
This project is licensed under the MIT License.
