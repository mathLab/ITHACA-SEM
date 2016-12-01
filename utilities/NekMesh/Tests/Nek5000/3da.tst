<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Nek5000 3D annulus</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list 3da.rea5000.gz 3da.xml:xml:test</parameters>
    <files>
        <file description="Input File">3da.rea5000.gz</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
