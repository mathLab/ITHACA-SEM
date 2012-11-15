<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Convergence: Taylor Vortex IMEXOrder2 dt=0.01</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>Test_TaylorVor_dt1.xml</parameters>
    <files>
        <file description="Session File">Test_TaylorVor_dt1.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">5.63502e-06</value>
            <value variable="v" tolerance="1e-12">5.81694e-06</value>
            <value variable="p" tolerance="1e-12">0.000232451</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">3.8882e-06</value>
            <value variable="v" tolerance="1e-12">4.91056e-06</value>
            <value variable="p" tolerance="1e-12">0.00016247</value>
        </metric>
    </metrics>
</test>