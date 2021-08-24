<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, pressure perturbation to test RiemannInvariant CBC supersonic, parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>--use-scotch Perturbation_M15_square_CBC_par.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File"> Perturbation_M15_square_CBC_par.xml</file>
        <file description="Restart File"> Perturbation_M15_square_CBC_par.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0910004</value>
            <value variable="rhou" tolerance="1e-12">40.244</value>
            <value variable="rhov" tolerance="1e-12">0.0038368</value>
            <value variable="E" tolerance="1e-12">23023.8</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.364797</value>
            <value variable="rhou" tolerance="1e-12">161.179</value>
            <value variable="rhov" tolerance="1e-12">0.431771</value>
            <value variable="E" tolerance="1e-12">92248.6</value>
        </metric>
    </metrics>
</test>


