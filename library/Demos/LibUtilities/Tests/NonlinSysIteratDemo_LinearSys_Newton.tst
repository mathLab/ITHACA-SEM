<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>NonlinSysIteratDemo for a small nonlinear problem of dimension 4</description>
    <executable>NonlinSysIteratDemo</executable>
    <parameters>NonlinSysIteratDemo_LinearSys_Newton.xml</parameters>
    <files>
        <file description="Session File">NonlinSysIteratDemo_LinearSys_Newton.xml</file>
    </files>

    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">9.99999965e-01</value>
            <value variable="v" tolerance="1e-12">9.26079326e-07</value>
        </metric>
    </metrics>
</test>
