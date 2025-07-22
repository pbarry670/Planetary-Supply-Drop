import numpy as np
import pandas as pd
from vpython import canvas, sphere, arrow, color, vector, distant_light, graph, gcurve, box, rate
import matplotlib.pyplot as plt

def runViz1():
    # Create a canvas
    scene1 = canvas(title="Simple VPython Animation", width=800, height=600)

    # Create a sphere
    ball = sphere(canvas=scene1, pos=vector(-5, 0, 0), radius=0.5, color=color.red)

    # Create a floor
    floor = box(canvas=scene1, pos=vector(0, -0.5, 0), size=vector(10, 0.1, 5), color=color.green)

    # Animation loop
    velocity = vector(0.1, 0, 0)  # Ball's velocity
    count = 1
    hasReachedCount = True
    while hasReachedCount:
        rate(50)  # Control the animation speed (50 frames per second)
        ball.pos += velocity  # Update the ball's position

        # Reverse direction if the ball hits the edge
        if ball.pos.x > 5 or ball.pos.x < -5:
            velocity.x = -velocity.x

        count = count + 1
        if (count >= 500):
            hasReachedCount = False

    print("Viz 1 complete.")
    

def runViz2():

    # Create a canvas
    scene2 = canvas(title="Simple VPython Animation", width=800, height=600)

    # Create a sphere
    ball = sphere(canvas=scene2, pos=vector(-5, 0, 0), radius=0.5, color=color.blue)

    # Create a floor
    floor = box(canvas=scene2, pos=vector(0, -0.5, 0), size=vector(10, 0.1, 5), color=color.green)

    # Animation loop
    velocity = vector(0.1, 0, 0)  # Ball's velocity
    count = 1
    hasReachedCount = True
    while hasReachedCount:
        rate(50)  # Control the animation speed (50 frames per second)
        ball.pos += velocity  # Update the ball's position

        # Reverse direction if the ball hits the edge
        if ball.pos.x > 5 or ball.pos.x < -5:
            velocity.x = -velocity.x

        count = count + 1
        if (count >= 500):
            hasReachedCount = False

    print("Viz 2 complete.")         

    

if __name__ == "__main__":

    print("Beginning.")
    runViz1()
    runViz2()
    print("End.")

