<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, pressure perturbation to test RiemannInvariant CBC, parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>--use-scotch Perturbation_M05_square_CBC_par.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File"> Perturbation_M05_square_CBC_par.xml</file>
        <file description="Restart File"> Perturbation_M05_square_CBC_par.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0910004</value>
            <value variable="rhou" tolerance="1e-12">13.4146</value>
            <value variable="rhov" tolerance="1e-12">0.00194307</value>
            <value variable="E" tolerance="1e-12">15113.8</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.364899</value>
            <value variable="rhou" tolerance="1e-12">53.726</value>
            <value variable="rhov" tolerance="1e-12">0.0856712</value>
            <value variable="E" tolerance="1e-12">60509.4</value>
        </metric>
    </metrics>
</test>


