<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, pressure perturbation to test RiemannInvariant CBC (back wave), parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>--use-scotch Perturbation_M05_circle_CBC_back_par.xml</parameters>
    <processes>6</processes>
    <files>
        <file description="Session File"> Perturbation_M05_circle_CBC_back_par.xml</file>
        <file description="Restart File"> Perturbation_M05_circle_CBC_back_par.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.080562</value>
            <value variable="rhou" tolerance="1e-12">11.8759</value>
            <value variable="rhov" tolerance="1e-12">0.00074691</value>
            <value variable="E" tolerance="1e-12">13380.1</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.364844</value>
            <value variable="rhou" tolerance="1e-12">53.7008</value>
            <value variable="rhov" tolerance="1e-12">0.0737579</value>
            <value variable="E" tolerance="1e-12">60483.6</value>
        </metric>
    </metrics>
</test>


