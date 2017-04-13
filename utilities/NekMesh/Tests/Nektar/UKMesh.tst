<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Nektar triangular mesh using older-style rea</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list UKMesh.rea.gz:rea:scalar UKMesh.xml:xml:test</parameters>
    <files>
        <file description="Input File">UKMesh.rea.gz</file>
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
