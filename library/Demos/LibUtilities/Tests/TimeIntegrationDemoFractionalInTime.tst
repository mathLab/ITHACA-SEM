<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test for time integration schemes</description>
    <executable>TimeIntegrationDemo</executable>
    <parameters>--dof 100 --timesteps 100 --method FractionalInTime --order 1 --test</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-8">8.80111e-08</value>
        </metric>
    </metrics>
</test>
