
import gym
from gym import spaces
import numpy as np
import simulator  # Import the C++ module

class BloodSuctionEnv(gym.Env):
    def __init__(self, grid_size=10, max_steps=200):
        super(BloodSuctionEnv, self).__init__()
        self.grid_size = grid_size
        self.max_steps = max_steps
        self.sim = simulator.Simulator(grid_size)

        self.action_space = spaces.Discrete(4)  # 4 actions: up, down, left, right
        self.observation_space = spaces.Box(low=0, high=1, shape=(grid_size, grid_size), dtype=np.float32)

        self.reset()

    def reset(self):
        self.sim.reset()
        self.current_step = 0
        return self.sim.get_state()

    def step(self, action):
        self.sim.step(action)
        state = self.sim.get_state()
        reward = self._calculate_reward(state)
        done = (np.sum(state) == 0) or (self.current_step >= self.max_steps)
        self.current_step += 1
        return state, reward, done, {}

    def _calculate_reward(self, state):
        return np.sum(1 - state)  # Reward for each cleared particle

    def render(self, mode='human'):
        import matplotlib.pyplot as plt
        grid = self.sim.get_state()
        plt.imshow(grid, cmap='Reds', interpolation='nearest')
        plt.show()
