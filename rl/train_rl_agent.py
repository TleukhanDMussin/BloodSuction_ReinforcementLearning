
import gym
from stable_baselines3 import PPO
from env.blood_suction_env import BloodSuctionEnv

# Create the environment
env = BloodSuctionEnv(grid_size=10, max_steps=200)

# Instantiate the PPO agent
model = PPO("MlpPolicy", env, verbose=1)

# Train the agent
model.learn(total_timesteps=100000)

# Save the trained model
model.save("ppo_blood_suction")

# Test the trained model
obs = env.reset()
done = False
while not done:
    action, _states = model.predict(obs, deterministic=True)
    obs, reward, done, _ = env.step(action)
    env.render()
