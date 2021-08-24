<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, pressure perturbation to test RiemannInvariant CBC supersonic (back wave), parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>--use-scotch Perturbation_M15_circle_CBC_back_par.xml</parameters>
    <processes>6</processes>
    <files>
        <file description="Session File"> Perturbation_M15_circle_CBC_back_par.xml</file>
        <file description="Restart File"> Perturbation_M15_circle_CBC_back_par.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0805624</value>
            <value variable="rhou" tolerance="1e-12">35.6279</value>
            <value variable="rhov" tolerance="1e-12">0.000129061</value>
            <value variable="E" tolerance="1e-12">20382.9</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.36433</value>
            <value variable="rhou" tolerance="1e-12">161.024</value>
            <value variable="rhov" tolerance="1e-12">0.0142433</value>
            <value variable="E" tolerance="1e-12">92155.9</value>
        </metric>
    </metrics>
</test>


