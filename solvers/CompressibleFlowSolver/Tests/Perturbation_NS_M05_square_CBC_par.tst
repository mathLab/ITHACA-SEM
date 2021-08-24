<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, pressure perturbation to test RiemannInvariant CBC, Navier-Stokes equations parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>--use-scotch Perturbation_NS_M05_square_CBC_par.xml</parameters>
    <processes>8</processes>
    <files>
        <file description="Session File"> Perturbation_NS_M05_square_CBC_par.xml</file>
        <file description="Restart File"> Perturbation_NS_M05_square_CBC_par.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">4.72522e-06</value>
            <value variable="rhou" tolerance="1e-12">0.000826153</value>
            <value variable="rhov" tolerance="1e-12">0.000972476</value>
            <value variable="E" tolerance="1e-12">9463.8</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.000219959</value>
            <value variable="rhou" tolerance="1e-12">0.0544127</value>
            <value variable="rhov" tolerance="1e-12">0.0262103</value>
            <value variable="E" tolerance="1e-12">37883.9</value>
        </metric>
    </metrics>
</test>


