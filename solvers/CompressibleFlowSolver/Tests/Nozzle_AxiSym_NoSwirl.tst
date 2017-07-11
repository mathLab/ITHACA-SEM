<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, axi-symmetric nozzle without swirl</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Nozzle_AxiSym_NoSwirl.xml</parameters>
    <files>
        <file description="Session File">Nozzle_AxiSym_NoSwirl.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">3.08375</value>
            <value variable="rhou" tolerance="1e-12">1.81859</value>
            <value variable="rhov" tolerance="1e-12">85.4653</value>
            <value variable="E" tolerance="1e-12">623440</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.26773</value>
            <value variable="rhou" tolerance="1e-12">1.812</value>
            <value variable="rhov" tolerance="1e-12">47.6855</value>
            <value variable="E" tolerance="1e-12">262359</value>
        </metric>
    </metrics>
</test>


