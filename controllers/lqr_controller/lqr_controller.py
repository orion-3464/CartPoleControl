from controller import Robot
from controller import Supervisor
import numpy as np
import math

m = 0.1
M = 0.135
l = 0.25
g = 9.81

K = [-10.0000, -8.3608, -32.6688, -5.4278]
    
k_e = 6.0    
k_p = 0.8    
k_d = 0.8    

MAX_FORCE = 30.0

ALPHA_LQR = 0.30
ALPHA_SWINGUP= 0.55

THETA_CAPTURE = 0.45 
THETA_ABANDON = 1.2

DURATION = 5.0
times = []

def run_robot():
    angles = []
    
    robot=Supervisor()
    dt=robot.getBasicTimeStep() / 1000.0
    timestep_ms=int(robot.getBasicTimeStep())


    cart_sensor=robot.getDevice("cart_sensor")
    cart_sensor.enable(timestep_ms)

    pole_sensor = robot.getDevice("pole_sensor")
    pole_sensor.enable(timestep_ms)

    cart_motor = robot.getDevice("cart_motor")
    cart_motor.setPosition(float('inf'))
    cart_motor.setForce(0.0)

    x_prev = 0.0
    theta_raw_prev = 0.0
    v_filtered = 0.0
    omega_filtered = 0.0
    lqr_active = False
    first_step = True

    while robot.step(timestep_ms) != -1:
        t = robot.getTime()
        if t >=DURATION:
            break
        x = cart_sensor.getValue()

        theta_raw = pole_sensor.getValue()
        theta = math.remainder(theta_raw + math.pi, 2 * math.pi)
        angles.append(theta)
        omega_raw = (theta_raw-theta_raw_prev)/dt

        v_raw = (x-x_prev)/dt
        
        if first_step:
            omega_filtered = omega_raw
            v_filtered = v_raw
            first_step = False
        
        if not lqr_active:
            alpha = ALPHA_SWINGUP
            if abs(theta) < THETA_CAPTURE:
                lqr_active = True
                omega_filtered = omega_raw
                v_filtered = v_raw
        else:
            alpha = ALPHA_LQR
            if abs(theta) > THETA_ABANDON:
                lqr_active = False
                omega_filtered = omega_raw
                v_filtered = v_raw
       

        omega_filtered = alpha*omega_filtered + (1.0-alpha)*omega_raw
        v_filtered = alpha * v_filtered + (1.0-alpha) * v_raw

        x_prev = x
        theta_raw_prev = theta_raw

        if lqr_active:
            u = -(K[0]*x+K[1]*v_filtered+K[2]*theta+ K[3]*omega_filtered)

        else:
            E = (2.0/3.0)*m*l**2*omega_filtered**2 + 1.5*m*g*l*math.cos(theta)
            E_star = 1.5*m*g*l
            E_error = E-E_star

            u_swing_accel =  k_e*E_error*omega_filtered*math.cos(theta)
            u_center_accel = -k_p*x-k_d*v_filtered
            x_ddot_des =  u_swing_accel+u_center_accel

            c  = math.cos(theta)
            s  = math.sin(theta)
            mass_term = (M + 2*m - (27.0/16.0)*m*c**2) * x_ddot_des
            gravity_term = (27.0/16.0)*m*g*s*c
            centrifugal_term = -1.5*m*l*omega_filtered**2*s

            u = mass_term + gravity_term + centrifugal_term

        u_clamped = max(min(u, MAX_FORCE), -MAX_FORCE)
        cart_motor.setForce(u_clamped)
        
    robot.simulationSetMode(Supervisor.SIMULATION_MODE_PAUSE)
    #robot.simulationSetMode(Supervisor.SIMULATION_MODE_RESET)
    angles = np.array(angles)
    np.save("angles.npy", angles)

if __name__ == "__main__":
    run_robot()