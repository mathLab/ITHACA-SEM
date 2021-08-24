<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 100 --timesteps 100 --method IMEX --variant DIRK --order 3 --parameter 3 4</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.00429031</value>
        </metric>
    </metrics>
</test>
