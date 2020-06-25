import probtools



students    = 20000
instructors =  2500
classes     =  3750
departments =   120
meeting_schedules = [[1,3],[1,3],[0,2,4],[0,2,4],[0,2]]

class_cohorts = 8

in_class_base_rate = 0.0145
in_dept_base_rate =  0.116
in_frnd_base_rate_hi = 0.18
in_frnd_base_rate_low = in_frnd_base_rate_hi * 0.25
in_dept_broad_base_rate = 0.00158


#_poisson_computer = probtools.Poisson()
#typical_R0 = 4.0
#typical_contacts_per_day = 10.2
#incubation_period = lambda : _poisson_computer.draw(5.2)
#poisson = _poisson_computer.draw
