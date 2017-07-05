<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, axi-symmetric nozzle with swirl</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Nozzle_AxiSym_Swirl.xml</parameters>
    <files>
        <file description="Session File">Nozzle_AxiSym_Swirl.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">3.08186</value>
            <value variable="rhou" tolerance="1e-12">1.82356</value>
            <value variable="rhov" tolerance="1e-12">85.3026</value>
            <value variable="rhow" tolerance="1e-12">10.8479</value>
            <value variable="E" tolerance="1e-12">622940</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.2665</value>
            <value variable="rhou" tolerance="1e-12">1.80257</value>
            <value variable="rhov" tolerance="1e-12">47.6879</value>
            <value variable="rhow" tolerance="1e-12">31.9047</value>
            <value variable="E" tolerance="1e-12">262000</value>
        </metric>
    </metrics>
</test>


