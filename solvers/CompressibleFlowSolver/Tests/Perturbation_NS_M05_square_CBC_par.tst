<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, pressure perturbation to test RiemannInvariant CBC, Navier-Stokes equations parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Perturbation_NS_M05_square_CBC_par.xml</parameters>
    <processes>8</processes>
    <files>
        <file description="Session File"> Perturbation_NS_M05_square_CBC_par.xml</file>
        <file description="Restart File"> Perturbation_NS_M05_square_CBC_par.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">4.71508e-06</value>
            <value variable="rhou" tolerance="1e-12">0.00082326</value>
            <value variable="rhov" tolerance="1e-12">0.000969901</value>
            <value variable="E" tolerance="1e-12">9463.8</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.000219754</value>
            <value variable="rhou" tolerance="1e-12">0.054418</value>
            <value variable="rhov" tolerance="1e-12">0.0261735</value>
            <value variable="E" tolerance="1e-12">37876.8</value>
        </metric>
    </metrics>
</test>


