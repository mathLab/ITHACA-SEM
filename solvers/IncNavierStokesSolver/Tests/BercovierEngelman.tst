<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Bercovier-Engelman unsteady Stokes flow</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>BercovierEngelman.xml</parameters>
    <files>
        <file description="Session File">BercovierEngelman.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">1.5619e-16</value>
            <value variable="v" tolerance="1e-12">1.94683e-16</value>
            <value variable="p" tolerance="5e-12">3.7773e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1.55431e-15</value>
            <value variable="v" tolerance="1e-12">1.4988e-15</value>
            <value variable="p" tolerance="1e-12">7.92052e-14</value>
        </metric>
    </metrics>
</test>
