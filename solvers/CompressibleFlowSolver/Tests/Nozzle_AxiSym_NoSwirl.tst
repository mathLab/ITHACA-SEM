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
            <value variable="rho" tolerance="1e-12">3.05647</value>
            <value variable="rhou" tolerance="1e-12">2.3371</value>
            <value variable="rhov" tolerance="1e-12">102.54</value>
            <value variable="E" tolerance="1e-12">616238</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.26181</value>
            <value variable="rhou" tolerance="1e-12">2.82123</value>
            <value variable="rhov" tolerance="1e-12">60.5981</value>
            <value variable="E" tolerance="1e-12">260653</value>
        </metric>
    </metrics>
</test>


