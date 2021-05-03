import numpy as np
import datetime
import csv
import math
#_________________________________________________________________________________________
# Программа вычисляет эфемериду малой планеты (α,δ) на 6 моментов времени:
# 01.06, 02.06, 03.06, 04.06, 05.06, 06.06 2021 года
#_________________________________________________________________________________________

def main():
    a, e, inc, ω, Ω, M0, date_array, elem_date, sun_coords = data()

    time_array = np.zeros(shape=(len(date_array)))
    for i in range(len(date_array)):
        time_array[i] = get_numDays(date_array[i], elem_date)

    inc, ω, Ω, M0 =  map(get_radian, [inc, ω, Ω, M0])

    α_array = np.zeros(shape=(len(date_array)))
    δ_array = np.zeros(shape=(len(date_array)))

    for i in range(len(date_array)):
        E = get_solve_kepler(a, e, M0, time_array[i])
        θ = get_true_anomaly(e, E)
        r = get_distance(E, a, e)
        ξ, η = get_orbit_coords(r, θ)

        x, y, z = get_eclip_coords(θ, inc, ω, Ω, r)
        x_eq, y_eq, z_eq = get_equat_coords(x, y, z)

        x_geo, y_geo, z_geo, ρ_mod = get_geo_coords(x_eq,  y_eq,  z_eq,
                                                    sun_coords[i,0],
                                                    sun_coords[i,1],
                                                    sun_coords[i,2])

        α, δ = get_efemerids(x_geo, y_geo, z_geo, ρ_mod)
        α_array[i] = get_degrees(α)
        δ_array[i] = get_degrees(δ)

    with open('results.csv', 'w', newline='') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(['Date', 'α', 'δ'])
        for i in range(time_array.size):
            writer.writerow([time_array[i],
                             get_hours(α_array[i]),
                             get_dms(δ_array[i])])

#_________________________________________________________________________________________
# ИСХОДНЫЕ ДАННЫЕ ЗАДАЧИ:
#_________________________________________________________________________________________
def data():
    a   = 2.620400  # большая полуось
    e   = 0.033909  # экстренситет
    inc = 5.71856   # наклон к эклиптике
    ω   = 181.76897 # аргумент перицентра
    Ω   = 70.25117  # долгота восходящего узла
    M0  = 187.59678 # средняя аномалия на 08.04
    #                            X           Y           Z
    sun_coords = np.array([[0.33939741, 0.87668949, 0.38003789],
                           [0.32342179, 0.88191245, 0.38230192],
                           [0.30735294, 0.88688572, 0.38445787],
                           [0.29119539, 0.89160770, 0.38650504],
                           [0.27495378, 0.89607685, 0.38844276],
                           [0.25863282, 0.90029174, 0.39027038]])

    date_array = ['1 Jun 21', '2 Jun 21', '3 Jun 21', '4 Jun 21', '5 Jun 21', '6 Jun 21']
    elem_date  = '8 Apr 21'

    return(a, e, inc, ω, Ω, M0,
           date_array, elem_date,
           sun_coords)

def get_numDays(input_date, elem_date):
    dt1 = datetime.datetime.strptime(input_date, "%d %b %y")
    dt2 = datetime.datetime.strptime(elem_date,  "%d %b %y")
    return (dt1 - dt2).days
#_________________________________________________________________________________________

def get_radian(varDegrees):
    in_radian = np.pi / 180
    var  =  varDegrees * in_radian
    return var
#_________________________________________________________________________________________
# 1.ОПРЕДЕЛЕНИЕ ОРБИТАЛЬНЫХ КООРДИНАТ:
#_________________________________________________________________________________________
def get_solve_kepler(a, e, M0, Δt):
    ϰ = 0.017202
    n = ϰ * a**(-1.5)
    M = M0 + n * Δt

    E = M
    i = 1
    while i < 10:
        E = M + e * np.sin(E)
        i += 1
    return E


def get_true_anomaly(e, E):
    θ = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E/2))
    if θ < 0:
        θ += 2 * np.pi
    return θ


def get_distance(E, a, e):
    r = a * (1 - e * np.cos(E))
    return r


def get_orbit_coords(r, θ):
    ξ = r * np.cos(θ)
    η = r * np.sin(θ)
    return ξ, η

#_________________________________________________________________________________________
# 2. ВЫЧИСЛЕНИЕ ГЕЛИОЦЕНТРИЧЕСКИХ КООРДИНАТ:
#_________________________________________________________________________________________
def get_eclip_coords(θ, i, ω, Ω, r):
    u = θ + ω
    x = r * (np.cos(u) * np.cos(Ω) - np.sin(u) * np.sin(Ω) * np.cos(i))
    y = r * (np.cos(u) * np.sin(Ω) + np.sin(u) * np.cos(Ω) * np.cos(i))
    z = r * np.sin(u)  * np.sin(i)
    return x, y, z


def get_equat_coords(x, y, z):
    eclip = get_radian(23.436742)
    x_eq = x
    y_eq = y * np.cos(eclip) - z * np.sin(eclip)
    z_eq = y * np.sin(eclip) + z * np.cos(eclip)
    return x_eq, y_eq, z_eq

#_________________________________________________________________________________________
# 3. ВЫЧИСЛЕНИЕ ГЕОЦЕНТРИЧЕСКИХ КООРДИНАТ:
#_________________________________________________________________________________________
def get_geo_coords(x_eq,  y_eq,  z_eq,
                   x_sun, y_sun, z_sun):
    x_geo = x_eq + x_sun
    y_geo = y_eq + y_sun
    z_geo = z_eq + z_sun
    ρ_mod = (x_geo**2 + y_geo**2 + z_geo**2)**0.5
    return x_geo, y_geo, z_geo, ρ_mod

#_________________________________________________________________________________________
# 4. ВЫЧИСЛЕНИЕ (α,δ):
#_________________________________________________________________________________________
def get_efemerids(x_geo, y_geo, z_geo, ρ_mod):
    α = np.arctan(y_geo / x_geo)
    δ = np.arcsin(z_geo / ρ_mod)
    if α < 0:
        α += np.pi
    return α, δ


def get_degrees(varRadian):
    in_degrees = 180 / np.pi
    var = varRadian * in_degrees
    return var


def get_dms(varDegrees):
    m, s = divmod(varDegrees*3600, 60)
    d, m = divmod(m, 60)
    d, m, s = map(int, [d, m, s])
    return f'{d}° {m}m {s}s'


def get_hours(varDegrees):
    m, s = divmod(varDegrees*240, 60)
    h, m = divmod(m, 60)
    h, m, s = map(int, [h, m, s])
    return f'{h}h {m}m {s}s'
#_________________________________________________________________________________________
if __name__ == '__main__':
    main()
