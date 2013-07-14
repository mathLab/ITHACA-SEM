<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, pressure perturbation to test RiemannInvariant CBC (back wave), parallel</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Perturbation_M05_square_CBC_back_par.xml</parameters>
    <processes>3</processes>
    <files>
        <file description="Session File"> Perturbation_M05_square_CBC_back_par.xml</file>
        <file description="Restart File"> Perturbation_M05_square_CBC_back_par.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.091</value>
            <value variable="rhou" tolerance="1e-12">13.4146</value>
            <value variable="rhov" tolerance="1e-12">0.000682517</value>
            <value variable="E" tolerance="1e-12">15113.7</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.365414</value>
            <value variable="rhou" tolerance="1e-12">53.7064</value>
            <value variable="rhov" tolerance="1e-12">0.108518</value>
            <value variable="E" tolerance="1e-12">60531.6</value>
        </metric>
    </metrics>
</test>


