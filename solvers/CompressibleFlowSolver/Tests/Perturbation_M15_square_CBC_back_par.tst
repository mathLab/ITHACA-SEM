<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, pressure perturbation to test RiemannInvariant CBC supersonic (back wave), parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>--use-scotch Perturbation_M15_square_CBC_back_par.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File"> Perturbation_M15_square_CBC_back_par.xml</file>
        <file description="Restart File"> Perturbation_M15_square_CBC_back_par.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0910005</value>
            <value variable="rhou" tolerance="1e-12">40.244</value>
            <value variable="rhov" tolerance="1e-12">0.000296026</value>
            <value variable="E" tolerance="1e-12">23023.8</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.364383</value>
            <value variable="rhou" tolerance="1e-12">161.034</value>
            <value variable="rhov" tolerance="1e-12">0.0306006</value>
            <value variable="E" tolerance="1e-12">92165.5</value>
        </metric>
    </metrics>
</test>


